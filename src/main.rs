use std::f64::consts::PI;

use crate::fitting::{
    calculate_statistics, find_parameters_gradient_descent, find_parameters_mcmc,
};
use crate::halo::{Halo, m200_c200_to_rs_rhos};
use crate::hydrostatics::{
    CM_IN_KPC, GG, KM_IN_KPC, isothermal_abg_background, isothermal_core_collapse_background,
};
use crate::logging::{load_file, save_output};
use crate::plotting::{create_corner_plot, plot_function};

mod fitting;
mod halo;
mod hydrostatics;
mod logging;
mod plotting;

const S_IN_GYR: f64 = 3.154e16;
const G_IN_MSUN: f64 = 1.988e33;

fn main() {
    let data = vec![
        (0.15649748494408605, 15595734576138530000.0),
        (0.18722388772639217, 15290963190523656000.0),
        (0.2252243432550796, 14924613361774703000.0),
        (0.26936457583283363, 14424784086219026000.0),
        (0.32301280500582225, 13894659794863262000.0),
        (0.38799031821006913, 13091285860285037000.0),
        (0.4647330501352986, 12151824423121148000.0),
        (0.5571917155289376, 11103690653134041000.0),
        (0.6692414165028057, 9918955303017603000.0),
        (0.801419168108796, 8606124456897346000.0),
        (0.9595907722097119, 7014426686842890000.0),
        (1.1500508076747187, 5250966932011177000.0),
        (1.3771319397378698, 3511010519320841000.0),
        (1.6463993031339876, 2293590617222256000.0),
        (1.983419488324015, 1615272370298475800.0),
    ];

    // Found by grad descent:
    let sidm_fit_params = [
        6282772.676310997,
        2.8897601910357062,
        0.18438405302532834,
        54574.17044567845,
    ];

    let cdm_fit_params = [
        7585724.648997071,
        2.195234319751577,
        0.0,
        169116.20672504138,
    ];

    let mcmc: bool = true;
    let premade: Option<String> = Some(String::from("32_x_10k.mcmc"));
    let params: [f64; 4];

    if mcmc {
        let (m200_params, chain, likelihoods): ([f64; 4], Vec<[f64; 4]>, Vec<f64>);
        let steps = 10000;
        let burn_in = 1000;
        if let Some(filename) = premade {
            let mcmc_output = load_file(filename).unwrap();
            m200_params = mcmc_output.best_params;
            chain = mcmc_output.chain;
            likelihoods = mcmc_output.likelihoods;
        } else {
            (m200_params, chain, likelihoods) = find_parameters_mcmc(
                &data,
                //[1e8, 5.0, 0.5,false//bad guess
                [2e9, 10.0, 0.2, 5.5e4],
                None,
                steps,
                burn_in,
                32,
            )
            .unwrap();

            save_output(
                format!("32_x_{}k.mcmc", steps / 1000),
                m200_params,
                chain.clone(),
                likelihoods,
            )
            .unwrap();
        }

        calculate_statistics(&chain, &m200_params);

        let bounds = [[1e7, 1e11], [0.0, 20.0], [0.0, 1.0], [1e2, 1e6]];

        check_chain_behavior(&chain);

        create_corner_plot(
            &chain,
            &["m200_0", "c200_0", "tau", "rho_c"],
            "corner_plot.png",
            (burn_in as f64) / (steps as f64),
            &bounds,
        )
        .unwrap();

        params = {
            let mut params = m200_params.clone();
            let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);
            params[0] = rho_s;
            params[1] = r_s;
            params
        }
    } else {
        //let params = find_parameters_gradient_descent(&data, [2e7, 4.0, 0.2, 5e5], None);
        params = sidm_fit_params;
    }

    let t = 10.0;
    let t_c = t / params[2];
    // [G rho] = km^2 kpc M_sun^-1 s^-2 * M_sun kpc^-3 = km^2 kpc^-2 s^-2
    // sigma_m = gyr^-1 M_sun^-1 kpc^2 kpc s km^-1 = (s/gyr) kpc^2 (kpc/km) M_sun^-1
    let mut sigma_m = 150.0
        / (0.75 * t_c * params[0] * params[1] * (4.0 * PI * GG * params[0]).sqrt())
        * (KM_IN_KPC / S_IN_GYR); // kpc^2 / M_sun
    sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // cm^2 / g

    println!("sigma_m = {sigma_m} \nt_c = {t_c}");

    let fit = isothermal_core_collapse_background(
        1e4,
        params[0],
        params[1],
        params[2],
        (0.1 * data[0].0, 10.0 * data.last().unwrap().0),
        params[3],
        false,
    );

    let legend_text = format!(
        "rho_s_0: {:.4e}\nr_s_0: {:.4e}\ntau: {}\nrho_c: {:.4e}",
        params[0], params[1], params[2], params[3]
    );
    println!("{}", &legend_text);

    plot_function(
        &fit.0,
        &fit.1,
        "fit.png",
        "Fit To Observed Profile",
        "r (arcmin)",
        "n_H (num / cm^2)",
        Some(legend_text),
        Some(&data),
    )
    .unwrap();

    let halo = Halo::NFW(params[0], params[1]);
    let r200 = halo.r200().unwrap();
    let m200 = halo.m200().unwrap();
    println!(
        "r200 = {}, m200 = {:.4e}, c200 = {} \ndeviation = {}",
        halo.r200().unwrap(),
        halo.m200().unwrap(),
        halo.c200().unwrap(),
        halo.deviation().unwrap()
    );
}

fn check_chain_behavior(chain: &[[f64; 4]]) {
    for i in 0..4 {
        let values: Vec<f64> = chain.iter().map(|p| p[i]).collect();
        let mean: f64 = values.iter().sum::<f64>() / values.len() as f64;
        let variance: f64 =
            values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
        let std = variance.sqrt();

        // Compare std to mean
        let relative_std = std / mean.abs();
        println!(
            "Param[{}]: mean = {:.6e}, std = {:.6e}, std/mean = {:.6e}",
            i, mean, std, relative_std
        );

        // For a narrow peak, std/mean should be small but > 0
        if relative_std < 1e-6 {
            println!("Warning: Chain may be stuck, not mixing!");
        } else if relative_std < 0.01 {
            println!("Info: Very narrow posterior (well-constrained parameter)");
        }
    }
}

fn plot_distribution() {
    let target_init_num_col_density = 1e19;
    let target_init_col_density: f64 =
        target_init_num_col_density * hydrostatics::M_PROTON * hydrostatics::CM_IN_KPC.powi(2);
    let temperature = 1e4;
    let mut sound_speed_squared = temperature / (hydrostatics::MP_OVER_KB); // kpc^2 / s^2
    sound_speed_squared *= hydrostatics::KM_IN_KPC.powi(2); // km^2 / s^2

    // r_s = c_s / sqrt(4pi G rho_c)
    // uniform sphere approx:
    // Sigma ~ 2*r_s*rho_c = 2 c_s sqrt(rho_c / 4pi G) => rho_c ~ 4pi G * (Sigma / 2c_s)^2 = pi G Sigma^2 / c_s^2
    // rho_c (kpc^-3) ~ pi G * (n_sigma (cm^-2) * MP * CMinKPC^2)^2 / c_s^2
    // Units: [rho_c] = M_sun kpc^-3 = [G * sigma^2 / c_s^2] = km^2 kpc M_sun^-1 s^-2 * M_sun^2 kpc^-4 * [c_s]^-2 = km^2 kpc^-3 M_sun s^-2 * [c_s]^-2
    // => [c_s^2] = km^2 kpc^-3 M_sun s^-2 / M_sun kpc^-3 = km^2 / s^2

    let rho_center_approx =
        PI * hydrostatics::GG * (target_init_col_density).powi(2) / sound_speed_squared;
    let scale_radius =
        sound_speed_squared.sqrt() / (4.0 * PI * hydrostatics::GG * rho_center_approx);
    println!("Center rho set to: {rho_center_approx}");

    // isothermal_abg_background(
    //     temperature,
    //     1.0,
    //     3.0,
    //     1.0,
    //     0.0e7,
    //     3.0,
    //     (1e-1, 10.0*scale_radius),
    //     rho_center_approx,
    // );

    isothermal_core_collapse_background(
        temperature,
        3e7,
        3.0,
        0.2,
        (1e-1, 10.0 * scale_radius),
        rho_center_approx,
        true,
    );
}
