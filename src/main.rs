use corner_plot::*;
use ensemble_mcmc::*;
use plotters::prelude::*;
use rayon::prelude::*;
use std::env;
use std::f64::consts::PI;
use std::fs;
use std::path::Path;
use std::path::PathBuf;

use crate::constants::*;
use crate::contour::get_3d_contour;
use crate::fitting::{Cloud9MCMCCore, Data};
use crate::halo::t_cross_section;
use crate::halo::{
    Halo, deviation, init_diemer_joyce, m200_c200_to_rs_rhos, mass_concentration_relation,
    rs_rhos_to_m200_c200,
};
use crate::hydrostatics::parametic_core_collapse;
use crate::hydrostatics::relhic_temperature_and_slope;
use crate::hydrostatics::{
    abg_background, core_collapse_background, evolution_profile, instability_profile, is_stable,
    relhic_neutral_fraction, relhic_temperature,
};
use crate::plot_utils::build::*;
use crate::plotting::*;
use crate::utils::*;

mod concentration_table;
mod constants;
mod contour;
// mod corner_plot;
mod fitting;
mod halo;
mod hydrostatics;
mod logging;
mod plot_utils;
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

    let args: Vec<usize> = env::args()
        .skip(1)
        .map(|a| a.parse::<usize>().expect("Expect a non-negative integer"))
        .collect();
    dbg!(&args);

    let mcmc_plots: bool = true;
    let prior = fitting::Prior::None;
    //let prior = fitting::Prior::MassConcentrationRelation(halo::McrSource::DiemerJoyce2019);

    let burn_in = 1000;
    let num_walkers = args[0]; //512;
    let steps = args[1]; //100000;
    let settings = MCMCSettings {
        burn_in,
        num_walkers,
        num_steps: steps,
        ..Default::default()
    };

    dbg!(settings);

    let fixed_cross_section = match args.len() == 2 {
        true => None,
        false => Some(args[2] as f64),
    };
    dbg!(fixed_cross_section);

    let file_name = make_file_name(num_walkers, steps, &prior, &fixed_cross_section);
    let data_path = String::from("data/") + &file_name;

    // run with let bounds = [[1e1, 5e9], [0.0, 9.0], [0.0, 1.0], [8e4, 1.6e5]];
    let bounds = [[1e9, 5e9], [0.0, 7.0], [0.0, 1.0], [8e4, 1.6e5]];
    let log_scales = vec![false, false, false, false]; // [m200, c200, tau, rho_c]

    //let font: FontDesc<'static> = ("sans-serif", 12).into_font(); //Paper
    let font: FontDesc<'static> = ("Times New Roman", 25).into_font(); //Presentations

    // likelihood_slice_profile(
    //     &data,
    //     &data_y_err_bar,
    //     [1e10, 6.5, 0.35, 3.5e4],
    //     0,
    //     [1e8, 1e30],
    //     fitting::Prior::MassConcentrationRelation(halo::McrSource::DuttonMaccio2014),
    // );

    // Active Code

    // create_cross_section_deviation_relation_plot(
    //     num_walkers,
    //     steps,
    //     &prior,
    //     &[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0],
    //     ("sans-serif", 30).into_font(),
    // )
    // .unwrap();
    // svg_to_pdf("figures/cross_section_vs_deviation_stable", 25.0).unwrap();

    let params: Vec<f64>;
    if mcmc_plots {
        let output: MCMCOutput;
        if Path::new(&(data_path.clone() + ".mcmc")).exists() {
            output = MCMCOutput::load(&(data_path.clone() + ".mcmc")).unwrap();
        } else {
            println!("Running MCMC!");
            let mcmc_core = Cloud9MCMCCore::init(data.clone(), prior, bounds, fixed_cross_section);
            output = mcmc(&mcmc_core, settings);

            println!("Gelman-Rubin R-hat statistics: {:?}", output.gelman_rubin);

            let converged = output.gelman_rubin.iter().all(|&r| r < 1.1);

            if !converged {
                println!("Warning: Chains may not have converged (R-hat > 1.1)");
                println!("Consider running more steps or tuning step sizes.");
            } else {
                println!("Chains appear to have converged (R-hat < 1.1)");
            }

            output.save(&(data_path.clone() + ".mcmc")).unwrap();
            output.save_as_json(&(data_path.clone() + ".json")).unwrap();
        }

        if !Path::new(&(data_path.clone() + "_stable.mcmc")).exists()
            || !Path::new(&(data_path.clone() + "_unstable.mcmc")).exists()
        {
            filter_stable(output.clone(), data_path.clone(), fixed_cross_section);
        }

        let stable_output = MCMCOutput::load(&(data_path.clone() + "_stable_bulk.mcmc")).unwrap();
        // let stable_chain = stable_output.chain.clone();
        // let unstable_output = MCMCOutput::load(&(data_path.clone() + "_unstable.mcmc")).unwrap();

        // filter_low_mass_island(unstable_output, (data_path.clone() + "_unstable"));

        // cdm_vs_sidm_fit_plot(&data);

        // instability_plot();

        // let num_stable = stable_output.chain.len();
        // let num_unstable = unstable_output.chain.len();

        // get_3d_contour(&stable_output.chain, &bounds, "cells_stable");
        // // get_3d_contour(&unstable_output.chain, &bounds, "cells_unstable");

        // filter_low_mass_island(stable_output.clone(), data_path.clone() + "_stable");

        // let island_output = MCMCOutput::load(&(data_path.clone() + "_stable_island.mcmc")).unwrap();
        // let num_island = island_output.chain.len();

        // dbg!(num_stable - num_island);
        // dbg!(num_unstable + num_island);
        // dbg!((num_unstable + num_island) as f64 / (num_stable + num_unstable) as f64);

        // let rhosrs_repack = island_output.chain[0..5]
        //     .iter()
        //     .map(|p: &Vec<f64>| {
        //         let (rs, rhos) = m200_c200_to_rs_rhos(p[0], p[1]);
        //         vec![rhos, rs, p[2], p[3]]
        //     })
        //     .collect();
        // instability_showcase(
        //     &rhosrs_repack,
        //     (1e1, 1e8),
        //     &relhic_temperature_and_slope,
        //     &relhic_neutral_fraction,
        //     ("sans-serif", 35).into_font(),
        // );

        // filter_cdm_island(stable_output.clone(), data_path.clone() + "_stable");

        // dbg!(cdmisland_output.chain.len());

        // create_mcr_deviation_plot(
        //     &stable_output.chain,
        //     &(String::from("figures/deviation_plot_") + &file_name + "_stable"),
        //     &[-8.0, 1.0],
        //     &[],
        //     font.clone(),
        // )
        // .unwrap();
        // create_mcr_deviation_plot(
        //     &unstable_output.chain,
        //     &(String::from("figures/deviation_plot_") + &file_name + "_unstable"),
        //     &[-8.0, 1.0],
        //     &[],
        //     font.clone(),
        // )
        // .unwrap();

        let corner_plot_format = CornerPlotFormat {
            axis_font: ("Times New Roman", 25).into_font(),
            legend_font: ("Times New Roman", 20).into_font(),
            log_scales: Some(log_scales),
            hist_bins: 75,
            contour_bins: 75,
            x_label_height: 60,
            y_label_width: 180,
            margin: 50,
            legend_widths: Some(vec![380, 220, 250, 350]),
            tex: true,
            do_hist_legend: false,
            ..Default::default()
        };
        create_corner_plot(
            &stable_output.chain,
            &[],
            // &["M₂₀₀", "c₂₀₀", "ρ꜀"],
            &["$M_{200}$", "$c_{200}$", "$\\tau$", "$\\rho_c$"],
            &(PathBuf::from(format!("figures/corner_plot_{file_name}.svg"))),
            &bounds,
            // &[bounds[0], bounds[1], bounds[3]],
            corner_plot_format.clone(),
        )
        .unwrap();
        svg_to_png(&format!("figures/corner_plot_{file_name}"), 25.0, true).unwrap();
        // create_corner_plot(
        //     &unstable_output.chain,
        //     &[],
        //     &["M₂₀₀", "c₂₀₀", "τ", "ρ꜀"],
        //     &(PathBuf::from(format!("figures/corner_plot_{file_name}_unstable.svg"))),
        //     &bounds,
        //     corner_plot_format,
        // )
        // .unwrap();
        // svg_to_png(&format!("figures/corner_plot_{file_name}_unstable")).unwrap();

        // check_chain_behavior(&chain);

        // get_3d_contour(&chain, &bounds);

        // let mcmc_sigma = {
        //     let params = output.best_params.clone();
        //     let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);

        //     let mut t_sigma_m = 150.0 * params[2]
        //         / (0.75 * r_s * rho_s * (4.0 * PI * GG * rho_s).sqrt())
        //         * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
        //     t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

        //     t_sigma_m
        // };

        // let (grad_descent_fit, grad_descent_sigma, grad_descent_deviation) = {
        //     let mut params = sidm_fit_params.clone();

        //     let mut t_sigma_m = 150.0 * params[2]
        //         / (0.75 * params[0] * params[1] * (4.0 * PI * GG * params[0]).sqrt())
        //         * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
        //     t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

        //     let (m200, c200) = rs_rhos_to_m200_c200(params[1], params[0]);
        //     params[0] = m200;
        //     params[1] = c200;

        //     let dev = deviation(m200, c200, halo::McrSource::DiemerJoyce2019);
        //     (params, t_sigma_m, dev)
        // };

        // let mut pruned_chain = Vec::new();
        // let mut marked_points = Vec::new();

        // for i in 0..chain.len() {
        //     let params = chain[i].clone();

        //     let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);
        //     let mut t_sigma_m = 150.0 * params[2]
        //         / (0.75 * rho_s * r_s * (4.0 * PI * GG * rho_s).sqrt())
        //         * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
        //     t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

        //     let sigma_m = t_sigma_m / 10.0;

        //     if sigma_m > 5e3 {
        //         marked_points.push(params.clone());
        //     }

        //     let dev = deviation(params[0], params[1], halo::McrSource::DiemerJoyce2019);
        //     if dev > -3.0 {
        //         pruned_chain.push(params)
        //     }
        // }
        // dbg!(marked_points.len());

        let cross_sec_chain: Vec<Vec<f64>> = output
            .chain
            .iter()
            .map(|params: &Vec<f64>| {
                let deviation = deviation(params[0], params[1], halo::McrSource::DiemerJoyce2019);
                let t_sigma_m = t_cross_section(params);
                // dbg!(t_sigma_m / 10.0, deviation);
                vec![t_sigma_m / 10.0, deviation]
            })
            .collect();

        let mini_corner_plot_format = CornerPlotFormat {
            axis_font: ("Times New Roman", 70).into_font(),
            legend_font: ("Times New Roman", 50).into_font(),
            log_scales: Some(vec![true, false]),
            x_label_height: 160,
            y_label_width: 370,
            margin: 80,
            hist_bins: 75,
            contour_bins: 75,
            legend_widths: Some(vec![250, 300]),
            tex: true,
            do_hist_legend: false,
            ..Default::default()
        };
        create_corner_plot(
            &cross_sec_chain,
            &[],
            &["$\\sigma/m$", "$\\textrm{Deviation}$"],
            &(PathBuf::from(format!("figures/mini_corner_plot_{file_name}.svg"))),
            &[[1e0, 5e4], [-7.0, 1.0]],
            mini_corner_plot_format,
        )
        .unwrap();
        svg_to_png(&format!("figures/mini_corner_plot_{file_name}"), 70.0, true).unwrap();

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
        //         let corner_plot_format = CornerPlotFormat {
        //             font: ("sans-serif", 35).into_font(),
        //             log_scales: Some(log_scales),
        //             hist_bins: 75,
        //             contour_bins: 75,
        //             x_label_height: 80,
        //             y_label_width: 140,
        //             ..Default::default()
        //         };
        //         create_corner_plot(
        //             &chain,
        //             &[],
        //             &["M₂₀₀", "C₂₀₀", "τ", "ρ꜀"],
        //             &(String::from("figures/corner_plot_") + &file_name),
        //             &bounds,
        //             corner_plot_format,
        //         )
        //         .unwrap();
        //     }
        //     Some(_) => {
        //         let corner_plot_format = CornerPlotFormat {
        //             font: font.clone(),
        //             log_scales: Some(vec![log_scales[0], log_scales[1], log_scales[3]]),
        //             ..Default::default()
        //         };
        //         create_corner_plot(
        //             &chain,
        //             &[],
        //             &["M₂₀₀", "C₂₀₀", "ρ꜀"],
        //             &(String::from("figures/corner_plot_") + &file_name),
        //             &[bounds[0], bounds[1], bounds[3]],
        //             corner_plot_format,
        //         )
        //         .unwrap();
        //     }
        // }

        // let m200_params = output.best_params.clone();
        let m200_params_idx = stable_output
            .log_likelihoods
            .iter()
            .enumerate()
            .max_by(|(_idxa, lla), (_idxb, llb)| lla.partial_cmp(llb).unwrap())
            .unwrap()
            .0;
        let m200_params = stable_output.chain[m200_params_idx].clone();

        params = match fixed_cross_section {
            None => {
                let mut params = m200_params;
                let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);
                params[0] = rho_s;
                params[1] = r_s;
                params
            }
            Some(cross_section) => {
                const AGE: f64 = 10.0;
                let (r_s, rho_s) = m200_c200_to_rs_rhos(m200_params[0], m200_params[1]);

                // cm^2 g^-1 Gyr kpc Ms kpc^-3 (km^2 kpc Ms^-1 s^-2 Ms kpc^-3)^0.5 = cm^2 g^-1 Gyr Ms kpc^-2 km kpc^-1 s^-1 = (cm/kpc)^2 (km/kpc) (Gyr/s) (Ms/g)
                let tau =
                    (0.75 * cross_section * AGE * r_s * rho_s * (4.0 * PI * GG * rho_s).sqrt()
                        / 150.0)
                        * S_IN_GYR
                        * G_IN_MSUN
                        / (KM_IN_KPC * CM_IN_KPC.powi(2));

                vec![rho_s, r_s, tau, m200_params[2]]
            }
        }
    } else {
        // params = find_parameters_gradient_descent(&data, sidm_fit_params, None, false);
        // params = find_parameters_gradient_descent(&data, [2e7, 4.0, 0.5, 5e5], None, false);
        // params = sidm_fit_params;
        panic!("Make sure you mean to do this");
    }

    dbg!(params.clone());

    // instability_profile(
    //     &params,
    //     (1e4, 1e7),
    //     &relhic_temperature,
    //     &relhic_neutral_fraction,
    //     font.clone(),
    // );
    // svg_to_pdf("figures/instability").unwrap();

    // evolution_profile(
    //     &params,
    //     &relhic_temperature,
    //     &relhic_neutral_fraction,
    //     &data,
    //     font.clone(),
    // );
    // svg_to_pdf("figures/evolution").unwrap();

    // let t = 10.0;
    // let t_c = t / params[2];
    // // [G rho] = km^2 kpc M_sun^-1 s^-2 * M_sun kpc^-3 = km^2 kpc^-2 s^-2
    // // sigma_m = gyr^-1 M_sun^-1 kpc^2 kpc s km^-1 = (s/gyr) kpc^2 (kpc/km) M_sun^-1
    // let mut sigma_m = 150.0
    //     / (0.75 * t_c * params[0] * params[1] * (4.0 * PI * GG * params[0]).sqrt())
    //     * (KM_IN_KPC / S_IN_GYR); // kpc^2 / M_sun
    // sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // cm^2 / g

    // println!("sigma_m = {sigma_m} \nt_c = {t_c}");

    // let halo = Halo::NFW(params[0], params[1]);

    // let fit = core_collapse_background(
    //     relhic_temperature,
    //     relhic_neutral_fraction,
    //     params[0],
    //     params[1],
    //     params[2],
    //     Some(1e20), //Some(params[3]),
    //     (1e-5 * INNER_BOUND, halo.r_crit()),
    //     true,
    // );

    // let halo = Halo::NFW(params[0], params[1]);
    // println!(
    //     "r200 = {}, m200 = {:.4e}, c200 = {} \ndeviation = {}",
    //     halo.r200().unwrap(),
    //     halo.m200().unwrap(),
    //     halo.c200().unwrap(),
    //     halo.deviation().unwrap()
    // );

    // let legend_text = format!(
    //     "rho_s_0: {:.2e}  \nr_s_0: {:.2}  \ntau: {:.2}",
    //     params[0], params[1], params[2]
    // );
    // println!("{}", &legend_text);

    // plot_functions(
    //     &fit.0,
    //     &vec![fit.1],
    //     &(String::from("figures/fit_") + &file_name + ".svg"),
    //     "Fit To Observed Profile",
    //     "r (arcmin)",
    //     "n_H (num / cm^2)",
    //     vec![Some(legend_text)],
    //     font,
    //     vec![false],
    //     Some(&data.points),
    //     Some(&data.y_err),
    // )
    // .unwrap();
}

fn filter_stable(output: MCMCOutput, data_path: String, fixed_cross_section: Option<f64>) {
    let chain = output.chain;
    let log_likelihoods = output.log_likelihoods;

    let (stable, unstable): (Vec<_>, Vec<_>) = chain
        .into_par_iter()
        .zip(log_likelihoods.into_par_iter())
        .partition(|(params, _)| {
            let (rs, rhos) = m200_c200_to_rs_rhos(params[0], params[1]);
            const AGE: f64 = 10.0;
            let (tau, rho_c) = match fixed_cross_section {
                None => (params[2], params[3]),
                Some(cross_section) => (
                    (0.75 * cross_section * AGE * rs * rhos * (4.0 * PI * GG * rhos).sqrt()
                        / 150.0)
                        * S_IN_GYR
                        * G_IN_MSUN
                        / (KM_IN_KPC * CM_IN_KPC.powi(2)),
                    params[2],
                ),
            };

            is_stable(
                &[rhos, rs, tau, rho_c],
                &relhic_temperature_and_slope,
                &relhic_neutral_fraction,
            )
        });

    let (stable_chain, stable_ll): (Vec<Vec<f64>>, Vec<f64>) = stable.into_iter().unzip();
    let (unstable_chain, unstable_ll): (Vec<Vec<f64>>, Vec<f64>) = unstable.into_iter().unzip();

    // Best params
    let stable_best_idx = stable_ll
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let unstable_best_idx = unstable_ll
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let stable_best_params = stable_chain[stable_best_idx].clone();
    let unstable_best_params = unstable_chain[unstable_best_idx].clone();

    let stable_output = MCMCOutput {
        chain: stable_chain,
        log_likelihoods: stable_ll,
        best_params: stable_best_params,
        gelman_rubin: output.gelman_rubin.clone(),
    };
    let unstable_output = MCMCOutput {
        chain: unstable_chain,
        log_likelihoods: unstable_ll,
        best_params: unstable_best_params,
        gelman_rubin: output.gelman_rubin.clone(),
    };

    stable_output
        .save(&(data_path.clone() + "_stable.mcmc"))
        .unwrap();
    unstable_output
        .save(&(data_path.clone() + "_unstable.mcmc"))
        .unwrap();
}

fn filter_low_mass_island(output: MCMCOutput, data_path: String) {
    let chain = output.chain;
    let log_likelihoods = output.log_likelihoods;

    let (island, bulk): (Vec<_>, Vec<_>) = chain
        .into_par_iter()
        .zip(log_likelihoods.into_par_iter())
        .partition(|(params, _)| params[0] < 1e9);

    let (island_chain, island_ll): (Vec<Vec<f64>>, Vec<f64>) = island.into_iter().unzip();
    let (bulk_chain, bulk_ll): (Vec<Vec<f64>>, Vec<f64>) = bulk.into_iter().unzip();

    // Best params
    let island_best_idx = island_ll
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let bulk_best_idx = bulk_ll
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let island_best_params = island_chain[island_best_idx].clone();
    let bulk_best_params = bulk_chain[bulk_best_idx].clone();

    let island_output = MCMCOutput {
        chain: island_chain,
        log_likelihoods: island_ll,
        best_params: island_best_params,
        gelman_rubin: output.gelman_rubin.clone(),
    };
    let bulk_output = MCMCOutput {
        chain: bulk_chain,
        log_likelihoods: bulk_ll,
        best_params: bulk_best_params,
        gelman_rubin: output.gelman_rubin.clone(),
    };

    island_output
        .save(&(data_path.clone() + "_island.mcmc"))
        .unwrap();
    bulk_output
        .save(&(data_path.clone() + "_bulk.mcmc"))
        .unwrap();
}

fn filter_cdm_island(output: MCMCOutput, data_path: String) {
    let chain = output.chain;
    let log_likelihoods = output.log_likelihoods;

    let (island, bulk): (Vec<_>, Vec<_>) = chain
        .into_par_iter()
        .zip(log_likelihoods.into_par_iter())
        .partition(|(params, _)| {
            t_cross_section(params) / 10.0 < 7.0
                && deviation(params[0], params[1], halo::McrSource::DiemerJoyce2019) < -5.0
        });

    let (island_chain, island_ll): (Vec<Vec<f64>>, Vec<f64>) = island.into_iter().unzip();
    let (bulk_chain, bulk_ll): (Vec<Vec<f64>>, Vec<f64>) = bulk.into_iter().unzip();

    // Best params
    let island_best_idx = island_ll
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let bulk_best_idx = bulk_ll
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let island_best_params = island_chain[island_best_idx].clone();
    let bulk_best_params = bulk_chain[bulk_best_idx].clone();

    let island_output = MCMCOutput {
        chain: island_chain,
        log_likelihoods: island_ll,
        best_params: island_best_params,
        gelman_rubin: output.gelman_rubin.clone(),
    };
    let bulk_output = MCMCOutput {
        chain: bulk_chain,
        log_likelihoods: bulk_ll,
        best_params: bulk_best_params,
        gelman_rubin: output.gelman_rubin.clone(),
    };

    island_output
        .save(&(data_path.clone() + "_cdmisland.mcmc"))
        .unwrap();
    bulk_output
        .save(&(data_path.clone() + "_sidmbulk.mcmc"))
        .unwrap();
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
