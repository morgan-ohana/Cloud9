use std::f64::consts::PI;
use std::fs;
use std::path::Path;

use crate::constants::*;
use crate::fitting::{
    calculate_gelman_rubin, calculate_statistics, find_parameters_gradient_descent,
    find_parameters_mcmc, likelihood_slice_profile, split_chains,
};
use crate::halo::{Halo, m200_c200_to_rs_rhos, rs_rhos_to_m200_c200};
use crate::hydrostatics::{isothermal_abg_background, isothermal_core_collapse_background};
use crate::logging::{load_file, save_output, save_output_json};
use crate::plotting::{create_chain_trace_plots, create_corner_plot, plot_function};

mod constants;
mod fitting;
mod halo;
mod hydrostatics;
mod logging;
mod plotting;

fn ensure_dir_exists(path: &str) -> Result<(), std::io::Error> {
    let path = Path::new(path);

    if !path.exists() {
        // Create directory and all parent directories if needed
        fs::create_dir_all(path)?;
    }

    Ok(())
}

fn main() {
    ensure_dir_exists("data").unwrap();
    ensure_dir_exists("trace_plots").unwrap();
    ensure_dir_exists("figures").unwrap();

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

    let data_y_err_bar = vec![
        (13732078146520720000.0, 17466755309837742000.0),
        (13427854263419597000.0, 17229242059185019000.0),
        (13003516161985374000.0, 16866179741318795000.0),
        (12461619020702566000.0, 16448095066290774000.0),
        (11695185329352206000.0, 15955195009096049000.0),
        (10716711258649842000.0, 15359790383829522000.0),
        (9581053603324199000.0, 14708083512691229000.0),
        (8205955659163100000.0, 13903084887830374000.0),
        (6793221053775176000.0, 12983171214351038000.0),
        (5337519460514788000.0, 11805617040307792000.0),
        (3880856839277999600.0, 10151147099502905000.0),
        (2370604088405835000.0, 8104671612570444000.0),
        (684370527200647700.0, 6345004385610885000.0),
        (4.4698143061e17, 5195521117287574000.0), // Lower bound for these two just y - 2(y_err_upper - y)
        (2.2656655652e17, 4312913465944557000.0), // ditto
    ];

    // Found by grad descent:
    let sidm_fit_params = [
        6282772.676310997,
        2.8897601910357062,
        0.18438405302532834,
        54574.17044567845,
    ];

    likelihood_slice_profile(
        &data,
        &data_y_err_bar,
        [1e10, 6.5, 0.35, 3.5e4],
        0,
        [1e8, 1e30],
    );

    let cdm_fit_params = [
        7585724.648997071,
        2.195234319751577,
        0.0,
        169116.20672504138,
    ];

    let mcmc: bool = true;

    let burn_in = 1000;
    let num_walkers = 512;
    let real_steps = 10000;
    let steps = real_steps + burn_in;

    let file_name = format!("{}_x_{}k", num_walkers, real_steps / 1000);
    let data_path = String::from("data/") + &file_name;

    let params: [f64; 4];
    let bounds = [[1e8, 1e16], [0.0, 20.0], [0.0, 1.0], [1e4, 1e6]];

    let premade: Option<String> = Some(data_path.clone() + ".mcmc");

    if mcmc {
        let (m200_params, chain, likelihoods): ([f64; 4], Vec<[f64; 4]>, Vec<f64>);
        if let Some(filename) = premade {
            let mcmc_output = load_file(filename).unwrap();
            m200_params = mcmc_output.best_params;
            chain = mcmc_output.chain;
            likelihoods = mcmc_output.likelihoods;
        } else {
            (m200_params, chain, likelihoods) = find_parameters_mcmc(
                &data,
                &data_y_err_bar,
                [2e9, 10.0, 0.2, 5.5e4],
                &bounds,
                steps,
                burn_in,
                num_walkers,
            );

            dbg!(chain.len());

            save_output(
                data_path.clone() + ".mcmc",
                m200_params,
                chain.clone(),
                likelihoods.clone(),
            )
            .unwrap();

            save_output_json(
                data_path.clone() + ".json",
                m200_params,
                chain.clone(),
                likelihoods.clone(),
            )
            .unwrap();
        }

        let mean_params = calculate_statistics(&chain, &m200_params);

        check_chain_behavior(&chain);

        dbg!(chain.len());
        dbg!(chain.len() as f64 / 32.0);

        // let chains = split_chains(&(m200_params, chain.clone(), likelihoods), real_steps);
        // create_chain_trace_plots(&chains).unwrap();

        let split_chains = split_chains(&(m200_params, chain.clone(), likelihoods), real_steps / 2);
        let r_hat = calculate_gelman_rubin(&split_chains);
        println!("Gelman-Rubin R-hat statistics: {:?}", r_hat);

        let converged = r_hat.iter().all(|&r| r < 1.1);

        if !converged {
            println!("Warning: Chains may not have converged (R-hat > 1.1)");
            println!("Consider running more steps or tuning step sizes.");
        } else {
            println!("Chains appear to have converged (R-hat < 1.1)");
        }

        let grad_descent_fit = {
            let mut params = sidm_fit_params.clone();
            let (m200, c200) = rs_rhos_to_m200_c200(params[1], params[0]);
            params[0] = m200;
            params[1] = c200;
            params
        };

        create_corner_plot(
            &chain,
            &[&grad_descent_fit, &mean_params, &m200_params],
            &["m200_0", "c200_0", "tau", "rho_c"],
            &(String::from("figures/corner_plot_") + &file_name),
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
        //params = find_parameters_gradient_descent(&data, [2e7, 4.0, 0.2, 5e5], None);
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
        UVB_TEMP,
        params[0],
        params[1],
        params[2],
        Some(params[3]),
        (0.1 * data[0].0, 10.0 * data.last().unwrap().0),
        false,
    );

    let halo = Halo::NFW(params[0], params[1]);
    println!(
        "r200 = {}, m200 = {:.4e}, c200 = {} \ndeviation = {}",
        halo.r200().unwrap(),
        halo.m200().unwrap(),
        halo.c200().unwrap(),
        halo.deviation().unwrap()
    );

    let legend_text = format!(
        "rho_s_0: {:.4e}\nr_s_0: {:.4e}\ntau: {}",
        params[0], params[1], params[2]
    );
    println!("{}", &legend_text);

    plot_function(
        &fit.0,
        &fit.1,
        &(String::from("figures/fit_") + &file_name + ".svg"),
        "Fit To Observed Profile",
        "r (arcmin)",
        "n_H (num / cm^2)",
        Some(legend_text),
        Some(&data),
    )
    .unwrap();
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
    let target_init_col_density: f64 = target_init_num_col_density * M_PROTON * CM_IN_KPC.powi(2);
    let temperature = UVB_TEMP;
    let mut sound_speed_squared = temperature / (MP_OVER_KB); // kpc^2 / s^2
    sound_speed_squared *= KM_IN_KPC.powi(2); // km^2 / s^2

    // r_s = c_s / sqrt(4pi G rho_c)
    // uniform sphere approx:
    // Sigma ~ 2*r_s*rho_c = 2 c_s sqrt(rho_c / 4pi G) => rho_c ~ 4pi G * (Sigma / 2c_s)^2 = pi G Sigma^2 / c_s^2
    // rho_c (kpc^-3) ~ pi G * (n_sigma (cm^-2) * MP * CMinKPC^2)^2 / c_s^2
    // Units: [rho_c] = M_sun kpc^-3 = [G * sigma^2 / c_s^2] = km^2 kpc M_sun^-1 s^-2 * M_sun^2 kpc^-4 * [c_s]^-2 = km^2 kpc^-3 M_sun s^-2 * [c_s]^-2
    // => [c_s^2] = km^2 kpc^-3 M_sun s^-2 / M_sun kpc^-3 = km^2 / s^2

    let rho_center_approx = PI * GG * (target_init_col_density).powi(2) / sound_speed_squared;
    let scale_radius = sound_speed_squared.sqrt() / (4.0 * PI * GG * rho_center_approx);
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

    let halo = Halo::NFW(3e7, 3.0);
    dbg!(halo.r_crit());
    dbg!(10.0 * scale_radius);
    let r_max = halo.r_crit();

    isothermal_core_collapse_background(
        temperature,
        3e7,
        3.0,
        0.2,
        Some(rho_center_approx),
        (1e-1, r_max),
        true,
    );
}
