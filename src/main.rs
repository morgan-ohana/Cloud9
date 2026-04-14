use plotters::prelude::*;
use std::f64::consts::PI;
use std::fs;
use std::path::Path;

use crate::constants::*;
use crate::contour::get_3d_contour;
use crate::corner_plot::create_corner_plot;
use crate::fitting::{Cloud9MCMCCore, Data};
use crate::halo::{
    Halo, deviation, init_diemer_joyce, m200_c200_to_rs_rhos, mass_concentration_relation,
    rs_rhos_to_m200_c200,
};
use crate::hydrostatics::{
    abg_background, core_collapse_background, evolution_profile, relhic_temperature,
};
use crate::logging::{load_file, save_output, save_output_json};
use crate::mcmc::{MCMCOutput, mcmc};
use crate::plotting::{
    create_cross_section_deviation_relation_plot, create_cross_section_plot,
    create_mcr_deviation_plot, plot_functions,
};
use crate::utils::*;

mod concentration_table;
mod constants;
mod contour;
mod corner_plot;
mod fitting;
mod halo;
mod hydrostatics;
mod logging;
mod mcmc;
mod plotting;
mod temperature;
mod utils;

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

    let _ = init_diemer_joyce("concentration_table.bin");

    let data = Data::init();

    // Control Switches

    let mcmc_plots: bool = true;
    let prior = fitting::Prior::None;
    //let prior = fitting::Prior::MassConcentrationRelation(halo::McrSource::DiemerJoyce2019);

    let burn_in = 1000;
    let num_walkers = 512;
    let steps = 100000;

    let fixed_cross_section = None; //Some(75.0);

    let file_name = make_file_name(num_walkers, steps, &prior, &fixed_cross_section);
    let data_path = String::from("data/") + &file_name;

    let bounds = [[1e8, 5e9], [0.0, 20.0], [0.0, 1.0], [8e4, 2e5]];
    let log_scales = [true, false, false, true]; // [m200, c200, tau, rho_c]

    //let font: FontDesc<'static> = ("sans-serif", 12).into_font(); //Paper
    let font: FontDesc<'static> = ("sans-serif", 25).into_font(); //Presentations

    let premade: Option<String> = Some(data_path.clone() + ".mcmc");

    // Cached Results
    // Found by grad descent:
    let sidm_fit_params_old = vec![
        6282772.676310997,
        2.8897601910357062,
        0.18438405302532834,
        54574.17044567845,
    ];

    let sidm_fit_params = vec![
        1699598.068854349,
        4.188800887645102,
        0.1721663180521211,
        116490.93147055767,
    ];

    // likelihood_slice_profile(
    //     &data,
    //     &data_y_err_bar,
    //     [1e10, 6.5, 0.35, 3.5e4],
    //     0,
    //     [1e8, 1e30],
    //     fitting::Prior::MassConcentrationRelation(halo::McrSource::DuttonMaccio2014),
    // );

    let cdm_fit_params = vec![
        7585724.648997071,
        2.195234319751577,
        0.0,
        169116.20672504138,
    ];

    let mut mcmc_fit_params = [1.9e9, 9.3, 1.5, 5.5e4]; //[2.95124e9, 8.925, 0.161084, 5.17612e4];

    {
        let (r_s, rho_s) = m200_c200_to_rs_rhos(mcmc_fit_params[0], mcmc_fit_params[1]);
        mcmc_fit_params[0] = rho_s;
        mcmc_fit_params[1] = r_s;
    }

    // Active Code

    cdm_vs_sidm_fit_plot(&data);

    // create_cross_section_deviation_relation_plot(
    //     num_walkers,
    //     steps,
    //     &prior,
    //     &[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0],
    //     font.clone(),
    // )
    // .unwrap();

    let params: Vec<f64>;
    if mcmc_plots {
        let output: MCMCOutput;
        if let Some(filename) = premade {
            output = load_file(filename).unwrap();
        } else {
            let mcmc_core = Cloud9MCMCCore::init(data.clone(), prior, bounds, fixed_cross_section);
            output = mcmc(&mcmc_core, steps, burn_in, num_walkers);

            save_output(data_path.clone() + ".mcmc", output.clone()).unwrap();
            save_output_json(data_path.clone() + ".json", output.clone()).unwrap();
        }

        println!("Gelman-Rubin R-hat statistics: {:?}", output.gelman_rubin);

        let converged = output.gelman_rubin.iter().all(|&r| r < 1.1);

        if !converged {
            println!("Warning: Chains may not have converged (R-hat > 1.1)");
            println!("Consider running more steps or tuning step sizes.");
        } else {
            println!("Chains appear to have converged (R-hat < 1.1)");
        }

        let chain = output.chain();

        check_chain_behavior(&chain);

        // get_3d_contour(&chain, &bounds);

        let mcmc_sigma = {
            let params = output.best_params.clone();
            let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);

            let mut t_sigma_m = 150.0 * params[2]
                / (0.75 * r_s * rho_s * (4.0 * PI * GG * rho_s).sqrt())
                * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
            t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

            t_sigma_m
        };
        dbg!(mcmc_sigma);

        let (grad_descent_fit, grad_descent_sigma, grad_descent_deviation) = {
            let mut params = sidm_fit_params.clone();

            let mut t_sigma_m = 150.0 * params[2]
                / (0.75 * params[0] * params[1] * (4.0 * PI * GG * params[0]).sqrt())
                * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
            t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

            let (m200, c200) = rs_rhos_to_m200_c200(params[1], params[0]);
            params[0] = m200;
            params[1] = c200;

            let dev = deviation(m200, c200, halo::McrSource::DiemerJoyce2019);
            (params, t_sigma_m, dev)
        };

        let mut pruned_chain = Vec::new();
        let mut marked_points = Vec::new();

        for i in 0..chain.len() {
            let params = chain[i].clone();

            let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);
            let mut t_sigma_m = 150.0 * params[2]
                / (0.75 * rho_s * r_s * (4.0 * PI * GG * rho_s).sqrt())
                * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
            t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

            let sigma_m = t_sigma_m / 10.0;

            if sigma_m > 5e3 {
                marked_points.push(params.clone());
            }

            let dev = deviation(params[0], params[1], halo::McrSource::DiemerJoyce2019);
            if dev > -3.0 {
                pruned_chain.push(params)
            }
        }
        dbg!(marked_points.len());

        // use rand::seq::SliceRandom;
        // let mut rng = rand::rng();
        // marked_points.shuffle(&mut rng);

        // let mut marked_points_arr: [&[f64; 4]; 100] = [&[0.0; 4]; 100];
        // let mut marked_dev: [f64; 100] = [0.0; 100];
        // for i in 0..100 {
        //     marked_points_arr[i] = &marked_points[i];
        //     marked_dev[i] = deviation(
        //         marked_points[i][0],
        //         marked_points[i][1],
        //         halo::McrSource::DiemerJoyce2019
        //     )
        // }

        // create_mcr_deviation_plot(
        //     &chain,
        //     &(String::from("figures/deviation_plot_") + &file_name),
        //     &[-8.0, 1.0],
        //     &[grad_descent_deviation],
        //     font.clone(),
        // )
        // .unwrap();

        // if let None = fixed_cross_section {
        //     create_cross_section_plot(
        //         &chain,
        //         &(String::from("figures/filtered_cross_section_plot_") + &file_name),
        //         &[1e-1, 1e7],
        //         &[grad_descent_sigma],
        //         font.clone(),
        //     )
        //     .unwrap();
        // }

        // match fixed_cross_section {
        //     None => {
        //         create_corner_plot(
        //             &chain,
        //             &[&grad_descent_fit],
        //             &["M₂₀₀", "C₂₀₀", "τ", "ρ꜀"],
        //             &log_scales,
        //             &(String::from("figures/corner_plot_") + &file_name),
        //             &bounds,
        //             font.clone(),
        //         )
        //         .unwrap();
        //     }
        //     Some(_) => {
        //         create_corner_plot(
        //             &chain,
        //             &[&grad_descent_fit],
        //             &["M₂₀₀", "C₂₀₀", "ρ꜀"],
        //             &[log_scales[0], log_scales[1], log_scales[3]],
        //             &(String::from("figures/corner_plot_") + &file_name),
        //             &[bounds[0], bounds[1], bounds[3]],
        //             font.clone(),
        //         )
        //         .unwrap();
        //     }
        // }

        params = match fixed_cross_section {
            None => {
                let mut params = output.best_params.clone(); // m200_params.clone();
                let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);
                params[0] = rho_s;
                params[1] = r_s;
                params
            }
            Some(cross_section) => {
                const AGE: f64 = 10.0;
                let (r_s, rho_s) =
                    m200_c200_to_rs_rhos(output.best_params[0], output.best_params[1]);

                // cm^2 g^-1 Gyr kpc Ms kpc^-3 (km^2 kpc Ms^-1 s^-2 Ms kpc^-3)^0.5 = cm^2 g^-1 Gyr Ms kpc^-2 km kpc^-1 s^-1 = (cm/kpc)^2 (km/kpc) (Gyr/s) (Ms/g)
                let tau =
                    (0.75 * cross_section * AGE * r_s * rho_s * (4.0 * PI * GG * rho_s).sqrt()
                        / 150.0)
                        * S_IN_GYR
                        * G_IN_MSUN
                        / (KM_IN_KPC * CM_IN_KPC.powi(2));

                vec![rho_s, r_s, tau, output.best_params[2]]
            }
        }
    } else {
        // params = find_parameters_gradient_descent(&data, sidm_fit_params, None, false);
        // params = find_parameters_gradient_descent(&data, [2e7, 4.0, 0.5, 5e5], None, false);
        params = sidm_fit_params;
    }

    dbg!(params.clone());

    evolution_profile(&params, &relhic_temperature, &data, font.clone());

    let t = 10.0;
    let t_c = t / params[2];
    // [G rho] = km^2 kpc M_sun^-1 s^-2 * M_sun kpc^-3 = km^2 kpc^-2 s^-2
    // sigma_m = gyr^-1 M_sun^-1 kpc^2 kpc s km^-1 = (s/gyr) kpc^2 (kpc/km) M_sun^-1
    let mut sigma_m = 150.0
        / (0.75 * t_c * params[0] * params[1] * (4.0 * PI * GG * params[0]).sqrt())
        * (KM_IN_KPC / S_IN_GYR); // kpc^2 / M_sun
    sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // cm^2 / g

    println!("sigma_m = {sigma_m} \nt_c = {t_c}");

    let halo = Halo::NFW(params[0], params[1]);

    let fit = core_collapse_background(
        relhic_temperature,
        params[0],
        params[1],
        params[2],
        Some(params[3]),
        (INNER_BOUND, halo.r_crit()),
        true,
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
        "rho_s_0: {:.2e}  \nr_s_0: {:.2}  \ntau: {:.2}",
        params[0], params[1], params[2]
    );
    println!("{}", &legend_text);

    plot_functions(
        &fit.0,
        &vec![fit.1],
        &(String::from("figures/fit_") + &file_name + ".svg"),
        "Fit To Observed Profile",
        "r (arcmin)",
        "n_H (num / cm^2)",
        vec![Some(legend_text)],
        font,
        vec![false],
        Some(&data.points),
        Some(&data.y_err),
    )
    .unwrap();
}

fn cdm_vs_sidm_fit_plot(data: &fitting::Data) {
    let sidm_output = load_file(String::from("data/512_x_100k.mcmc")).unwrap();
    let (sidm_rs, sidm_rhos) =
        m200_c200_to_rs_rhos(sidm_output.best_params[0], sidm_output.best_params[1]);
    let sidm_halo = Halo::NFW(sidm_rhos, sidm_rs);

    // let cdm_params = sidm_output
    //     .chain()
    //     .into_iter()
    //     .zip(sidm_output.log_likelihoods.clone())
    //     .filter(|(p, _)| p[2] < 0.05)
    //     .max_by(|(_, la), (_, lb)| la.partial_cmp(lb).unwrap())
    //     .map(|(p, _)| p)
    //     .unwrap();
    let cdm_params = load_file(String::from("data/512_x_100k_sigma=0.mcmc"))
        .unwrap()
        .best_params;
    let (cdm_rs, cdm_rhos) = m200_c200_to_rs_rhos(cdm_params[0], cdm_params[1]);
    let cdm_halo = Halo::NFW(cdm_rhos, cdm_rs);

    let collapse_params = sidm_output
        .chain()
        .into_iter()
        .zip(sidm_output.log_likelihoods)
        .filter(|(p, _)| p[2] > 0.95)
        .max_by(|(_, la), (_, lb)| la.partial_cmp(lb).unwrap())
        .map(|(p, _)| p)
        .unwrap();
    let (collapse_rs, collapse_rhos) = m200_c200_to_rs_rhos(collapse_params[0], collapse_params[1]);
    let collapse_halo = Halo::NFW(collapse_rhos, collapse_rs);

    let r_max = cdm_halo
        .r_crit()
        .max(sidm_halo.r_crit())
        .max(collapse_halo.r_crit()); //shared r_max so they share r_grid

    let cdm_fit = core_collapse_background(
        relhic_temperature,
        cdm_rhos,
        cdm_rs,
        0.0,
        Some(cdm_params[2]),
        (INNER_BOUND, r_max),
        false,
    );

    let tau = sidm_output.best_params[2];
    let sidm_fit = core_collapse_background(
        relhic_temperature,
        sidm_rhos,
        sidm_rs,
        tau,
        Some(sidm_output.best_params[3]),
        (INNER_BOUND, r_max),
        false,
    );

    let collapse_tau = collapse_params[2];
    let collapse_fit = core_collapse_background(
        relhic_temperature,
        collapse_rhos,
        collapse_rs,
        collapse_tau,
        Some(collapse_params[3]),
        (INNER_BOUND, r_max),
        false,
    );

    plot_functions(
        &cdm_fit.0,
        &vec![cdm_fit.1, sidm_fit.1, collapse_fit.1],
        "cdm_vs_sidm.svg",
        "Best CDM and SIDM fits",
        "r (arcmin)",
        "H₁ Column Density (cm⁻²)",
        vec![
            Some(String::from("CDM")),
            Some(format!("SIDM τ={tau:.2}")),
            Some(format!("SIDM τ={collapse_tau:.2}")),
        ],
        ("sans-serif", 25).into_font(),
        vec![false, true, true],
        Some(&data.points),
        Some(&data.y_err),
    )
    .unwrap()
}

fn check_chain_behavior(chain: &Vec<Vec<f64>>) {
    for i in 0..chain[0].len() {
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
