use corner_plot::{CornerPlotFormat, plot_histogram};
use ensemble_mcmc::MCMCOutput;
use plotters::coord::Shift;
use plotters::prelude::full_palette::ORANGE_700;
use plotters::prelude::*;
use rayon::prelude::*;

use crate::constants::*;
use crate::fitting;
use crate::halo::{Halo, McrSource, deviation, m200_c200_to_rs_rhos};
use crate::hydrostatics::core_collapse_background_3d_output;
use crate::hydrostatics::instability_showcase;
use crate::hydrostatics::is_stable;
use crate::hydrostatics::parametic_core_collapse;
use crate::hydrostatics::relhic_temperature_and_slope;
use crate::hydrostatics::{core_collapse_background, relhic_neutral_fraction, relhic_temperature};
// use crate::logging::load_file;
use crate::plot_utils::{build::*, frame::*, legend::*, utils::*};
use crate::utils::make_file_name;
use core::f64;
use std::f64::consts::PI;

/// 1.0e9 -> "1.0\times 10^{9}" (no $...$, so it composes into larger expressions)
pub fn sci_latex(x: f64) -> String {
    let exp = x.abs().log10().floor() as i32;
    let mantissa = x / 10f64.powi(exp);
    if (mantissa - 1.0).abs() < 1e-6 {
        format!("10^{{{exp}}}")
    } else {
        format!("{mantissa:.1}\\times 10^{{{exp}}}")
    }
}

pub fn fmt_num(num: &f64) -> String {
    if num.abs() <= 1e-100 {
        return "$0$".to_string();
    }
    if num.abs() >= 1000.0 || num.abs() <= 0.1 {
        format!("${}$", sci_latex(*num))
    } else {
        format!("${num:.1}$")
    }
}

pub fn plot_functions(
    x_points: &Vec<f64>,
    y_points: &Vec<Vec<f64>>,
    filename: &str,
    title: &str,
    xlabel: &str,
    ylabel: &str,
    legends: Vec<Option<String>>,
    font: FontDesc<'static>,
    dashed: Vec<bool>,
    data: Option<&Vec<(f64, f64)>>,
    data_y_err: Option<&Vec<(f64, f64)>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new(filename, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    draw_functions_on_area(
        &root, x_points, y_points, title, xlabel, ylabel, legends, font, dashed, data, data_y_err,
    )?;

    root.present()?;
    Ok(())
}

pub fn draw_functions_on_area<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    x_points: &Vec<f64>,
    y_points: &Vec<Vec<f64>>,
    title: &str,
    xlabel: &str,
    ylabel: &str,
    legends: Vec<Option<String>>,
    font: FontDesc<'static>,
    dashed: Vec<bool>,
    data: Option<&Vec<(f64, f64)>>,
    data_y_err: Option<&Vec<(f64, f64)>>,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;
    let x_min;
    let x_max;

    let data_spread = match data {
        Some(data_points) => Some(data_points.last().unwrap().0 - data_points[0].0),
        None => None,
    };

    match data {
        Some(data_points) => {
            x_min = 0.5 * data_points[0].0;
            x_max = 2.0 * data_points.last().unwrap().0;
        }
        // Some(data_points) => {
        //     // x_min = x_points[0];
        //     // x_max = *x_points.last().unwrap();
        //     x_min = data_points[0].0 - 0.1 * data_spread.unwrap();
        //     x_max = data_points.last().unwrap().0 + 0.1 * data_spread.unwrap();
        // }
        None => {
            x_min = x_points[0];
            x_max = *x_points.last().unwrap();
        }
    };

    let x_range = x_min..x_max;

    match data {
        Some(data_points) => {
            for i in 0..data_points.len() {
                if data_points[i].1 > y_max {
                    y_max = data_points[i].1
                }
                if data_points[i].1 < y_min {
                    y_min = data_points[i].1
                }
            }

            if let Some(y_err) = data_y_err {
                for i in 0..y_err.len() {
                    if y_err[i].1 > y_max {
                        y_max = y_err[i].1
                    }
                    if y_err[i].0 < 1e-10 {
                        continue;
                    }
                    if y_err[i].0 < y_min {
                        y_min = y_err[i].0
                    }
                }
            }

            for n in 0..y_points.len() {
                for i in 0..y_points[n].len() {
                    if !x_range.contains(&x_points[i]) {
                        continue;
                    }
                    if y_points[n][i] > y_max {
                        y_max = y_points[n][i]
                    }
                    if y_points[n][i] < y_min {
                        y_min = y_points[n][i]
                    }
                }
            }
        }
        None => {
            for n in 0..y_points.len() {
                for i in 0..y_points[n].len() {
                    if y_points[n][i] > y_max {
                        y_max = y_points[n][i]
                    }
                    if y_points[n][i] < y_min {
                        y_min = y_points[n][i]
                    }
                }
            }
        }
    }
    dbg!(y_min);

    // Buffer y_bounds
    let (y_min, y_max) = (
        (y_min + 1e-4)
            * match y_min.signum() {
                1.0 => 0.9,
                -1.0 => 1.1,
                _ => panic!("number has no sign, is probably NaN"),
            },
        y_max
            * match y_max.signum() {
                1.0 => 1.25,
                -1.0 => 0.9,
                _ => panic!("number has no sign, is probably NaN"),
            },
    );

    // println!("y_min = {:.3}, y_max={:.3}", y_min, y_max);

    dbg!(y_min);
    let y_range = y_min..y_max;

    let log_scale = false;

    let x_range = x_range.log_scale();
    let y_range = y_range.log_scale();
    let log_scale = true;

    let mut chart = ChartBuilder::on(area)
        .caption(title, ("Times New Roman", 40))
        .margin(30)
        .x_label_area_size(40) //30)
        .y_label_area_size(120) //60)
        .build_cartesian_2d(x_range, y_range)?;

    chart
        .configure_mesh()
        .x_desc(xlabel) // X-axis label
        .y_desc(ylabel) // Y-axis label
        .x_label_style(font.clone())
        .y_label_style(font.clone())
        .x_label_formatter(&fmt_num)
        .y_label_formatter(&fmt_num)
        .x_labels(3)
        .y_labels(15)
        .axis_style(BLACK.stroke_width(1))
        .draw()?;

    let mut plot_profiles: Vec<Vec<(f64, f64)>> = Vec::with_capacity(y_points.len());
    for n in 0..y_points.len() {
        plot_profiles.push(
            (0..x_points.len())
                .map(|i| (x_points[i], y_points[n][i]))
                .collect(),
        );
    }

    let dashed_count = dashed.iter().filter(|&&d| d).count().max(1);
    let mut dashed_index = 0;
    let mut legend_entries: Vec<LegendEntry<DB>> = Vec::new();
    for n in 0..y_points.len() {
        let color = if dashed[n] {
            let t = dashed_index as f64 / (dashed_count - 1).max(1) as f64;

            // Simple linear gradient: BLUE → RED
            let r = (255.0 * t) as u8;
            let g = 0;
            let b = (255.0 * (1.0 - t)) as u8;

            dashed_index += 1;

            RGBColor(r, g, b)
        } else {
            BLACK
        };

        let series = if dashed[n] {
            chart.draw_series(DashedLineSeries::new(
                plot_profiles[n].clone(),
                10,
                10,
                ShapeStyle {
                    color: color.mix(1.0),
                    filled: false,
                    stroke_width: 1,
                },
            ))?
        } else {
            chart.draw_series(LineSeries::new(plot_profiles[n].clone(), &color))?
        };

        if let Some(legend_text) = &legends[n] {
            legend_entries.push(LegendEntry {
                tex: legend_text.clone(),
                proxy: latex_to_proxy(legend_text),
                marker: if dashed[n] {
                    LegendMarker::Dashed(color, 5, 3)
                } else {
                    LegendMarker::Line(color)
                },
            });
        }
    }

    if let Some(data_points) = data {
        chart
            .draw_series(
                data_points
                    .iter()
                    .map(|point| Circle::new(*point, 5, &BLACK)),
            )
            .unwrap();
    }

    if let (Some(data_points), Some(y_err)) = (data, data_y_err) {
        assert_eq!(data_points.len(), y_err.len());

        let cap_width = 0.02; // fraction of x in log space (tweak this)

        chart.draw_series(
            data_points
                .iter()
                .zip(y_err.iter())
                .map(|(&(x, _y), &(yl, yh))| {
                    let mut elems = Vec::new();

                    // Vertical error bar
                    elems.push(PathElement::new(vec![(x, yl), (x, yh)], &RED));

                    let (x_left, x_right) = match log_scale {
                        true => (x * (1.0 - cap_width), x * (1.0 + cap_width)),
                        false => {
                            let cap_width = cap_width
                                * (data_points[data_points.len() - 1].0 - data_points[0].0);
                            (x - 0.5 * cap_width, x + 0.5 * cap_width)
                        }
                    };

                    elems.push(PathElement::new(vec![(x_left, yl), (x_right, yl)], &RED));
                    elems.push(PathElement::new(vec![(x_left, yh), (x_right, yh)], &RED));

                    elems
                })
                .flatten(),
        )?;
    }

    if !legend_entries.is_empty() {
        draw_legend(
            &chart.plotting_area().strip_coord_spec(),
            &legend_entries,
            ("sans-serif", 25).into_font(),
            LegendAnchor::UpperRight,
            10,
            None, // or Some(px) to force the box width
        )?;
    }

    let (x_ticks, y_ticks) = if log_scale {
        (log_ticks(x_min, x_max, 0), log_ticks(y_min, y_max, -1))
    } else {
        (linear_ticks(x_min, x_max, 3), linear_ticks(y_min, y_max, 5))
    };

    draw_ticks_top_and_right(
        &chart,
        &x_ticks,
        &y_ticks,
        (x_min, x_max),
        (y_min, y_max),
        5,
        BLACK.stroke_width(1),
    )?;

    //println!("Plot saved as {}", filename);
    Ok(())
}

pub fn instability_plot() {
    let stable = MCMCOutput::load("data/512_x_100k_stable.mcmc").unwrap();
    let unstable = MCMCOutput::load("data/512_x_100k_unstable.mcmc").unwrap();

    let mut stable_params = stable.best_params;
    let mut unstable_params = unstable.chain[1].clone();
    dbg!(&stable_params);
    dbg!(&unstable_params);

    (stable_params[1], stable_params[0]) = m200_c200_to_rs_rhos(stable_params[0], stable_params[1]);
    (unstable_params[1], unstable_params[0]) =
        m200_c200_to_rs_rhos(unstable_params[0], unstable_params[1]);

    dbg!(is_stable(
        &stable_params,
        &relhic_temperature_and_slope,
        &relhic_neutral_fraction
    ));
    dbg!(is_stable(
        &unstable_params,
        &relhic_temperature_and_slope,
        &relhic_neutral_fraction
    ));
    instability_showcase(
        &vec![stable_params, unstable_params],
        (1.5e3, 1e8),
        &relhic_temperature_and_slope,
        &relhic_neutral_fraction,
        ("Times New Roman", 30).into_font(),
    );
}

pub fn cdm_vs_sidm_fit_plot(data: &fitting::Data) {
    let mut params: Vec<Vec<f64>> = Vec::new();
    let mut labels = Vec::new();

    // SIDM BEST
    let sidm_output = MCMCOutput::load("data/512_x_100k_stable_bulk.mcmc").unwrap();
    // params.push(sidm_output.best_params.clone());
    // labels.push(Some(format!("SIDM τ={:.2}", sidm_output.best_params[2])));

    // CDM BEST
    let cdm_params = MCMCOutput::load("data/512_x_100k_sigma=0_stable_bulk.mcmc")
        .unwrap()
        .best_params;
    params.push(vec![cdm_params[0], cdm_params[1], 0.0, cdm_params[2]]);
    labels.push(Some(format!("$\\textrm{{NFW}}$")));

    // CORE FORMING
    let core_params = sidm_output
        .chain
        .iter()
        .zip(&sidm_output.log_likelihoods)
        .filter(|(p, _)| p[2] > 0.175 && p[2] < 0.225)
        .max_by(|(_, la), (_, lb)| la.partial_cmp(lb).unwrap())
        .map(|(p, _)| p)
        .unwrap();
    params.push(core_params.clone());
    labels.push(Some(format!(
        "$\\textrm{{SIDM}} \\, \\tau={:.2}$",
        core_params[2]
    )));

    // COLLAPSE
    let collapse_params = sidm_output
        .chain
        .iter()
        .zip(&sidm_output.log_likelihoods)
        .filter(|(p, _)| p[2] > 0.95)
        .max_by(|(_, la), (_, lb)| la.partial_cmp(lb).unwrap())
        .map(|(p, _)| p)
        .unwrap();
    params.push(collapse_params.clone());
    labels.push(Some(format!(
        "$\\textrm{{SIDM}} \\,\\tau={:.2}$",
        collapse_params[2]
    )));

    let mut r_max: f64 = 0.0;
    let mut rs_rhos: Vec<(f64, f64)> = Vec::new();
    for p in &params {
        let (m200, c200, tau) = (p[0], p[1], p[2]);
        let (rs, rhos) = m200_c200_to_rs_rhos(m200, c200);
        let halo = Halo::NFW(rs, rhos);
        let dev = deviation(m200, c200, McrSource::DiemerJoyce2019);
        let mut t_sigma_m = 150.0 * tau / (0.75 * rhos * rs * (4.0 * PI * GG * rhos).sqrt())
            * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
        t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

        r_max = r_max.max(halo.r_crit());

        rs_rhos.push((rs, rhos));
    }
    r_max *= 1e2;

    let (rs_min, rs_max) = rs_rhos
        .iter()
        .map(|(rs, _)| rs)
        .fold((f64::MAX, f64::MIN), |(min, max), &val| {
            (min.min(val), max.max(val))
        });

    let (rhos_min, rhos_max) = rs_rhos
        .iter()
        .map(|(_, rhos)| rhos)
        .fold((f64::MAX, f64::MIN), |(min, max), &val| {
            (min.min(val), max.max(val))
        });

    let mut gas_column_density = Vec::new();
    let mut ang_points = Vec::new();
    for i in 0..params.len() {
        let (tau, rhoc) = (params[i][2], params[i][3]);
        let (rs, rhos) = rs_rhos[i];

        let fit = core_collapse_background(
            relhic_temperature_and_slope,
            relhic_neutral_fraction,
            rhos,
            rs,
            tau,
            Some(rhoc),
            (INNER_BOUND, r_max),
            false,
        );
        if i == 0 {
            ang_points = fit.0;
        }
        gas_column_density.push(fit.1);
    }

    let root = SVGBackend::new("figures/cdm_vs_sidm_stable.svg", (1024, 768)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let (x_min, x_max) = (0.5 * data.points[0].0, 2.0 * data.points.last().unwrap().0);
    let x_range = x_min..x_max;

    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;
    for &(_, y) in &data.points {
        y_min = y_min.min(y);
        y_max = y_max.max(y);
    }
    for &(yl, yh) in &data.y_err {
        y_max = y_max.max(yh);
        if yl > 1e-10 {
            y_min = y_min.min(yl);
        }
    }
    for (_, col) in gas_column_density.iter().enumerate() {
        for (i, &y) in col.iter().enumerate() {
            if !x_range.contains(&ang_points[i]) {
                continue;
            }
            y_min = y_min.min(y);
            y_max = y_max.max(y);
        }
    }

    let (y_min, y_max) = (
        (y_min + 1e-4)
            * match y_min.signum() {
                1.0 => 0.9,
                -1.0 => 1.1,
                _ => panic!("number has no sign, is probably NaN"),
            },
        y_max
            * match y_max.signum() {
                1.0 => 1.25,
                -1.0 => 0.9,
                _ => panic!("number has no sign, is probably NaN"),
            },
    );

    let mut chart = ChartBuilder::on(&root)
        .caption("", ("Times New Roman", 40))
        .margin(30)
        .x_label_area_size(40)
        .y_label_area_size(100)
        .build_cartesian_2d((x_min..x_max).log_scale(), (y_min..y_max).log_scale())
        .unwrap();

    chart
        .configure_mesh()
        .x_desc("$r \\,(\\mathrm{arcmin})$")
        .y_desc("$\\mathrm{H}\\textsc{i}\\, \\textrm{Column Density}\\, (\\mathrm{cm}^{-2}$)")
        .x_label_style(("Times New Roman", 25).into_font())
        .y_label_style(("Times New Roman", 25).into_font())
        .x_label_formatter(&fmt_num)
        .y_label_formatter(&fmt_num)
        .x_labels(3)
        .y_labels(15)
        .axis_style(BLACK.stroke_width(3))
        .draw()
        .unwrap();

    let mut legend_entries: Vec<LegendEntry<SVGBackend>> = Vec::new();
    let colors = vec![BLACK, BLUE, RED];
    for (n, col) in gas_column_density.iter().enumerate() {
        let color = colors[n];
        let profile: Vec<(f64, f64)> = window_series(
            &ang_points
                .iter()
                .copied()
                .zip(col.iter().copied())
                .collect::<Vec<(f64, f64)>>(),
            (x_min, x_max),
            (y_min, y_max),
        );

        chart
            .draw_series(DashedLineSeries::new(
                profile,
                10,
                10,
                ShapeStyle {
                    color: color.mix(1.0),
                    filled: false,
                    stroke_width: 2,
                },
            ))
            .unwrap();

        if let Some(tex) = &labels[n] {
            legend_entries.push(LegendEntry {
                tex: tex.clone(),
                proxy: latex_to_proxy(tex) + "11", //manual buffer characters
                marker: LegendMarker::Dashed(color, 5, 3),
            });
        }
    }

    chart
        .draw_series(
            data.points
                .iter()
                .map(|p| Circle::new(*p, 5, Into::<ShapeStyle>::into(&ORANGE_700).filled())),
        )
        .unwrap();

    {
        let cap_width = 0.05; // fraction of x in log space
        chart
            .draw_series(data.points.iter().zip(data.y_err.iter()).flat_map(
                |(&(x, _), &(yl, yh))| {
                    let (x_left, x_right) = (x * (1.0 - cap_width), x * (1.0 + cap_width));
                    let bar_style = ORANGE_700.stroke_width(3);
                    let mut bar_series = vec![PathElement::new(vec![(x, yl), (x, yh)], bar_style)];

                    if yl > y_min {
                        bar_series.push(PathElement::new(
                            vec![(x_left, yl), (x_right, yl)],
                            bar_style,
                        ));
                    }
                    if yh < y_max {
                        bar_series.push(PathElement::new(
                            vec![(x_left, yh), (x_right, yh)],
                            bar_style,
                        ));
                    }
                    bar_series
                },
            ))
            .unwrap();
    }

    draw_legend(
        &chart.plotting_area().strip_coord_spec(),
        &legend_entries,
        ("Times New Roman", 25).into_font(),
        LegendAnchor::UpperRight,
        10,
        None,
    )
    .unwrap();

    let (x_ticks, y_ticks) = (log_ticks(x_min, x_max, 0), log_ticks(y_min, y_max, 0));
    draw_ticks_top_and_right(
        &chart,
        &x_ticks,
        &y_ticks,
        (x_min, x_max),
        (y_min, y_max),
        5,
        BLACK.stroke_width(3),
    )
    .unwrap();

    fit_comp_plot(data, params, labels);

    root.present().unwrap();

    svg_to_pdf("figures/cdm_vs_sidm_stable", 25.0).unwrap();
}

pub fn fit_comp_plot(data: &fitting::Data, params: Vec<Vec<f64>>, labels: Vec<Option<String>>) {
    let mut rs_rhos: Vec<(f64, f64)> = Vec::new();
    for p in &params {
        let (m200, c200, tau) = (p[0], p[1], p[2]);
        let (rs, rhos) = m200_c200_to_rs_rhos(m200, c200);
        let dev = deviation(m200, c200, McrSource::DiemerJoyce2019);
        let mut t_sigma_m = 150.0 * tau / (0.75 * rhos * rs * (4.0 * PI * GG * rhos).sqrt())
            * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
        t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

        rs_rhos.push((rs, rhos));
    }

    let (rs_min, rs_max) = rs_rhos
        .iter()
        .map(|(rs, _)| rs)
        .fold((f64::MAX, f64::MIN), |(min, max), &val| {
            (min.min(val), max.max(val))
        });

    let (rhos_min, rhos_max) = rs_rhos
        .iter()
        .map(|(_, rhos)| rhos)
        .fold((f64::MAX, f64::MIN), |(min, max), &val| {
            (min.min(val), max.max(val))
        });

    let (rmin, rmax) = (1e-2 * rs_min, 5e0 * rs_max);
    let (rhomin, rhomax) = (1e-2 * rhos_min, 5e2 * rhos_max);
    const R_GRID_NUM: usize = 1000;
    let inset_r_range: Vec<f64> = (0..R_GRID_NUM)
        .into_iter()
        .map(|i: usize| -> f64 {
            (rmin.ln() + (i as f64 / (R_GRID_NUM - 1) as f64) * (rmax.ln() - rmin.ln())).exp()
        })
        .collect();

    let mut gas_3d_density: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut gas_neutral_3d_density: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut dm_density = Vec::new();
    for i in 0..params.len() {
        let (tau, rhoc) = (params[i][2], params[i][3]);
        let (rs, rhos) = rs_rhos[i];

        let dm_density_func = parametic_core_collapse(rs, rhos, tau);
        let dm_rho_pts: Vec<(f64, f64)> = inset_r_range
            .iter()
            .map(|&r| (r, dm_density_func(r)))
            .collect();
        dm_density.push(window_series(&dm_rho_pts, (rmin, rmax), (rhomin, rhomax)));

        let gas_3d = core_collapse_background_3d_output(
            relhic_temperature_and_slope,
            relhic_neutral_fraction,
            rhos,
            rs,
            tau,
            Some(rhoc),
            &inset_r_range,
        );
        gas_3d_density.push(window_series(
            &inset_r_range
                .clone()
                .into_iter()
                .zip(gas_3d.1)
                .collect::<Vec<(f64, f64)>>(),
            (rmin, rmax),
            (rhomin, rhomax),
        ));
        gas_neutral_3d_density.push(window_series(
            &inset_r_range
                .clone()
                .into_iter()
                .zip(gas_3d.0)
                .collect::<Vec<(f64, f64)>>(),
            (rmin, rmax),
            (rhomin, rhomax),
        ));
    }

    let root = SVGBackend::new("figures/fit_comparison.svg", (1024, 768)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .caption("", ("Times New Roman", 40))
        .margin(30)
        .x_label_area_size(40)
        .y_label_area_size(100)
        .build_cartesian_2d((rmin..rmax).log_scale(), (rhomin..rhomax).log_scale())
        .unwrap();

    chart
        .configure_mesh()
        .x_desc("$r\\, (\\mathrm{kpc})$")
        .y_desc("$\\rho$ ($\\mathrm{M}_\\odot \\mathrm{kpc}^{-3}$)")
        .x_label_style(("Times New Roman", 25).into_font())
        .y_label_style(("Times New Roman", 25).into_font())
        .x_labels(3)
        .y_labels(3)
        .x_label_formatter(&fmt_num)
        .y_label_formatter(&fmt_num)
        .axis_style(BLACK.stroke_width(3))
        .draw()
        .unwrap();

    // Build/Plot density profiles

    let colors = vec![BLACK, BLUE, RED];
    let mut legend_entries: Vec<LegendEntry<SVGBackend>> = Vec::new();
    for i in 0..params.len() {
        let color = colors[i];
        chart
            .draw_series(LineSeries::new(
                dm_density[i].clone(),
                color.stroke_width(2),
            ))
            .unwrap();

        let tex = format!(
            "$\\textrm{{Dark Matter:}} \\, M_{{200}}$={}, $c_{{200}}$={}",
            fmt_num(&params[i][0]),
            fmt_num(&params[i][1]),
        );
        legend_entries.push(LegendEntry {
            proxy: latex_to_proxy(&tex),
            tex,
            marker: LegendMarker::Line(color),
        });
    }

    for i in 0..params.len() {
        let color = colors[i];

        chart
            .draw_series(DashedLineSeries::new(
                gas_3d_density[i].clone(),
                3,
                10,
                ShapeStyle {
                    color: color.mix(1.0),
                    filled: false,
                    stroke_width: 2,
                },
            ))
            .unwrap();
        let tex = format!("$\\textrm{{Gas:}}\\, \\rho_c$={}", fmt_num(&params[i][3]));
        legend_entries.push(LegendEntry {
            proxy: latex_to_proxy(&tex),
            tex,
            marker: LegendMarker::Dashed(color, 2, 4),
        });
    }
    for i in 0..params.len() {
        let color = colors[i];
        chart
            .draw_series(DashedLineSeries::new(
                gas_neutral_3d_density[i].clone(),
                10,
                10,
                ShapeStyle {
                    color: color.mix(1.0),
                    filled: false,
                    stroke_width: 2,
                },
            ))
            .unwrap();
    }

    let tex = format!("$\\textrm{{Neutral Gas}}$");
    let marker_colors = colors.clone();
    legend_entries.push(LegendEntry {
        proxy: latex_to_proxy(&tex),
        tex,
        marker: LegendMarker::Custom(Box::new(
            move |area: &DrawingArea<SVGBackend, Shift>,
                  y: i32,
                  x0: i32,
                  x1: i32|
                  -> Result<(), Box<dyn std::error::Error>> {
                const ON: i32 = 5;
                const OFF: i32 = 3;
                let mut x = x0;
                let mut k = 0usize;
                while x < x1 {
                    let xe = (x + ON).min(x1);
                    area.draw(&PathElement::new(
                        vec![(x, y), (xe, y)],
                        marker_colors[k % marker_colors.len()].stroke_width(2),
                    ))?;
                    k += 1;
                    x += ON + OFF;
                }
                Ok(())
            },
        )),
    });

    draw_legend(
        &chart.plotting_area().strip_coord_spec(),
        &legend_entries,
        ("Times New Roman", 25).into_font(),
        LegendAnchor::UpperRight,
        10,
        Some(420),
    )
    .unwrap();

    let (x_ticks, y_ticks) = (log_ticks(rmin, rmax, 0), log_ticks(rhomin, rhomax, 0));

    draw_ticks_top_and_right(
        &chart,
        &x_ticks,
        &y_ticks,
        (rmin, rmax),
        (rhomin, rhomax),
        5,
        BLACK.stroke_width(3),
    )
    .unwrap();

    root.present().unwrap();

    svg_to_pdf("figures/fit_comparison", 25.0).unwrap();
}

pub fn create_cross_section_deviation_relation_plot(
    num_walkers: usize,
    steps: usize,
    prior: &crate::fitting::Prior,
    cross_sections: &[f64],
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    let deviations: Vec<_> = cross_sections.to_vec().par_iter().map(|cross_section| {
        let file_name = make_file_name(num_walkers, steps, prior, &Some(*cross_section));

        if let Ok(output) = MCMCOutput::load(&(String::from("data/") + &file_name + "_stable_bulk.mcmc")) {
            println!("File loaded for sigma/m={cross_section}");
            let chain = output.chain;

            let mut mean_deviation = 0.0;
            let mut mean_square_deviation = 0.0;
            for params in &chain {
                let deviation = deviation(params[0], params[1], McrSource::DiemerJoyce2019);

                mean_deviation += deviation;
                mean_square_deviation += deviation.powi(2);
            }

            mean_deviation /= chain.len() as f64;
            mean_square_deviation /= chain.len() as f64;
            let var_deviation = mean_square_deviation - mean_deviation.powi(2);

            Some((*cross_section, mean_deviation, var_deviation))
        } else {
            println!(
                "Unable to find data for cross_section = {cross_section}! Has this MCMC been run?"
            );
            None
        }
    }).collect();

    let mut data: Vec<(f64, f64)> = Vec::with_capacity(deviations.len());
    let mut y_err: Vec<(f64, f64)> = Vec::with_capacity(deviations.len());
    for deviation_data in deviations.into_iter().flatten() {
        let (cross_section, mean_deviation, var_deviation): (f64, f64, f64) = deviation_data;
        data.push((cross_section, mean_deviation));
        y_err.push((
            mean_deviation - var_deviation.sqrt(),
            mean_deviation + var_deviation.sqrt(),
        ));
    }

    plot_functions(
        &Vec::new(),
        &Vec::new(),
        &(String::from("figures/cross_section_vs_deviation_stable_bulk.svg")),
        "Deviation Dependence on Cross Section",
        "Cross Section (cm² / g)",
        "Deviation",
        Vec::new(),
        font,
        Vec::new(),
        Some(&data),
        Some(&y_err),
    )?;

    Ok(())
}

pub fn create_mcr_deviation_plot(
    chain: &[Vec<f64>],
    output_path: &str,
    bounds: &[f64; 2],
    marked_values: &[f64],
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut deviation_chain = Vec::with_capacity(chain.len());
    let mut mean_deviation = 0.0;
    let mut mean_square_deviation = 0.0;
    for params in chain {
        let deviation = deviation(params[0], params[1], McrSource::DiemerJoyce2019);

        mean_deviation += deviation;
        mean_square_deviation += deviation.powi(2);
        deviation_chain.push(deviation);
    }

    mean_deviation /= chain.len() as f64;
    mean_square_deviation /= chain.len() as f64;
    let var_deviation = mean_square_deviation - mean_deviation.powi(2);

    dbg!(mean_deviation, var_deviation.sqrt());

    let output_path = String::from(output_path);
    let svg_path = &(output_path.clone() + ".svg");
    let pdf_path = &(output_path.clone() + ".pdf");
    // Create plot area
    let root = SVGBackend::new(svg_path, (1600, 1600)).into_drawing_area();
    root.fill(&WHITE)?;

    let corner_plot_format = CornerPlotFormat {
        font,
        ..Default::default()
    };

    plot_histogram(
        &root.margin(
            5,
            5 + corner_plot_format.x_label_height,
            5 + corner_plot_format.y_label_width,
            5,
        ),
        &deviation_chain,
        "deviation from cosmological median",
        (0, 0),
        1,
        &bounds,
        &marked_values,
        None,
        &corner_plot_format,
    )?;

    root.present()?;
    println!("Sigma plot saved to: {}.svg", output_path);

    let status = std::process::Command::new("inkscape")
        .args(&["--export-type=pdf", "--export-filename", pdf_path, svg_path])
        .status()?;

    if status.success() {
        Ok(())
    } else {
        Err(format!("Inkscape failed with status: {}", status).into())
    }
}

pub fn create_cross_section_plot(
    chain: &[Vec<f64>],
    output_path: &str,
    bounds: &[f64; 2],
    marked_values: &[f64],
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut t_sigma_m_chain = Vec::with_capacity(chain.len());
    let mut mean_log_sigma = 0.0;
    let mut mean_square_log_sigma = 0.0;
    for params in chain {
        let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);
        let mut t_sigma_m = 150.0 * params[2]
            / (0.75 * rho_s * r_s * (4.0 * PI * GG * rho_s).sqrt())
            * (KM_IN_KPC / S_IN_GYR); // Gyr kpc^2 / M_sun
        t_sigma_m *= CM_IN_KPC.powi(2) / G_IN_MSUN; // Gyr cm^2 / g

        mean_log_sigma += t_sigma_m.log10();
        mean_square_log_sigma += t_sigma_m.log10().powi(2);
        t_sigma_m_chain.push(t_sigma_m.log10());
    }

    mean_log_sigma /= chain.len() as f64;
    mean_square_log_sigma /= chain.len() as f64;
    let var_log_sigma = mean_square_log_sigma - mean_log_sigma.powi(2);

    dbg!(mean_log_sigma, var_log_sigma.sqrt());

    let bounds = [bounds[0].log10(), bounds[1].log10()];
    let marked_values = [marked_values[0].log10()];

    let output_path = String::from(output_path);
    let svg_path = &(output_path.clone() + ".svg");
    let pdf_path = &(output_path.clone() + ".pdf");
    // Create plot area
    let root = SVGBackend::new(svg_path, (1600, 1600)).into_drawing_area();
    root.fill(&WHITE)?;

    let corner_plot_format = CornerPlotFormat {
        font,
        ..Default::default()
    };

    plot_histogram(
        &root.margin(
            5,
            5 + corner_plot_format.x_label_height,
            5 + corner_plot_format.y_label_width,
            5,
        ),
        &t_sigma_m_chain,
        "log(t sigma/m)",
        (0, 0),
        1,
        &bounds,
        &marked_values,
        None, //Some(&[&prob_dens]),
        &corner_plot_format,
    )?;

    root.present()?;
    println!("Sigma plot saved to: {}.svg", output_path);

    let status = std::process::Command::new("inkscape")
        .args(&["--export-type=pdf", "--export-filename", pdf_path, svg_path])
        .status()?;

    if status.success() {
        Ok(())
    } else {
        Err(format!("Inkscape failed with status: {}", status).into())
    }
}
