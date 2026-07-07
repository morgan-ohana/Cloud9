use corner_plot::{CornerPlotFormat, plot_histogram};
use ensemble_mcmc::MCMCOutput;
use plotters::coord::ranged1d::ValueFormatter;
use plotters::prelude::*;
use rayon::prelude::*;

use crate::constants::*;
use crate::fitting;
use crate::halo::{Halo, McrSource, deviation, m200_c200_to_rs_rhos};
use crate::hydrostatics::instability_showcase;
use crate::hydrostatics::parametic_core_collapse;
use crate::hydrostatics::{core_collapse_background, relhic_neutral_fraction, relhic_temperature};
// use crate::logging::load_file;
use crate::utils::make_file_name;
use crate::utils::svg_to_pdf;
use std::f64::consts::PI;

fn fmt_num(num: &f64) -> String {
    if num.abs() <= 1e-100 {
        // Probably true 0
        return format!("0");
    }

    if num.abs() >= 1000.0 || num.abs() <= 0.1 {
        format!("{:.1e}", num)
    } else {
        format!("{:.1}", num)
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

fn draw_functions_on_area<DB: DrawingBackend>(
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
        //     // x_min = data_points[0].0 - 0.1 * data_spread.unwrap();
        //     // x_max = data_points.last().unwrap().0 + 0.1 * data_spread.unwrap();
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
                1.0 => 1.1,
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
        .caption(title, ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(60) //30)
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
        .y_labels(5)
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
    let mut make_legend = false;
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
            series.label(legend_text).legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle {
                        color: color.mix(1.0),
                        filled: false,
                        stroke_width: 2,
                    },
                )
            });

            make_legend = true
        }
    }
    if make_legend {
        // Configure and draw legend
        chart
            .configure_series_labels()
            .label_font(font)
            .position(SeriesLabelPosition::UpperRight)
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;
    }

    if let Some(data_points) = data {
        chart
            .draw_series(data_points.iter().map(|point| Circle::new(*point, 5, &RED)))
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

    let (x_ticks, y_ticks) = if log_scale {
        (log_ticks(x_min, x_max, 0), log_ticks(y_min, y_max, 0))
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

fn linear_ticks(min: f64, max: f64, n: usize) -> Vec<f64> {
    let step = (max - min) / (n as f64);
    // Round step to a "nice" number (1, 2, 2.5, 5 × power of 10)
    let magnitude = 10f64.powf(step.log10().floor());
    let nice_step = [1.0, 2.0, 2.5, 5.0, 10.0]
        .iter()
        .map(|&f| f * magnitude)
        .find(|&s| s >= step)
        .unwrap_or(magnitude);

    // Start at the first multiple of nice_step >= min
    let first = (min / nice_step).ceil() * nice_step;

    std::iter::successors(Some(first), |&v| {
        let next = v + nice_step;
        if next <= max + 1e-10 * nice_step {
            Some(next)
        } else {
            None
        }
    })
    .collect()
}

fn log_ticks(min: f64, max: f64, level: i32) -> Vec<f64> {
    assert!(min > 0.0);
    assert!(max > 0.0);
    let exp_min = min.log10().floor() as i32;
    let exp_max = max.log10().ceil() as i32;

    if level >= 0 {
        // Coarser than one-per-decade: tick every 10^level decades
        let stride = 10i32.pow(level as u32);
        let first_exp = (exp_min as f64 / stride as f64).ceil() as i32 * stride;
        std::iter::successors(Some(first_exp), |&e| Some(e + stride))
            .take_while(|&e| e <= exp_max)
            .map(|exp| 10f64.powi(exp))
            .filter(|&v| v >= min * (1.0 - 1e-10) && v <= max * (1.0 + 1e-10))
            .collect()
    } else {
        // Sub-decade: step = 10^(exp + level + 1) within each decade.
        // Decade boundaries are generated as k=0 of each decade — no double-counting.
        (exp_min..=exp_max)
            .flat_map(move |exp| {
                let decade_start = 10f64.powi(exp);
                let step = 10f64.powi(exp + level + 1);
                let next_decade = 10f64.powi(exp + 1);
                (0..)
                    .map(move |k| decade_start + k as f64 * step)
                    .take_while(move |&v| v < next_decade * (1.0 - 1e-10))
            })
            .filter(|&v| v >= min * (1.0 - 1e-10) && v <= max * (1.0 + 1e-10))
            .collect()
    }
}

fn draw_ticks_top_and_right<
    DB: DrawingBackend,
    X: Ranged<ValueType = f64> + ValueFormatter<f64>,
    Y: Ranged<ValueType = f64> + ValueFormatter<f64>,
>(
    chart: &ChartContext<DB, Cartesian2d<X, Y>>,
    x_ticks: &[f64],
    y_ticks: &[f64],
    (x_min, x_max): (f64, f64),
    (y_min, y_max): (f64, f64),
    tick_size: i32, // pixels; positive = inward
    style: ShapeStyle,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let coord_spec = chart.as_coord_spec();

    let (x_px_range, y_px_range) = (
        coord_spec.get_x_axis_pixel_range(),
        coord_spec.get_y_axis_pixel_range(),
    );

    let (x_offset, y_offset) = (
        x_px_range.start.min(x_px_range.end),
        y_px_range.start.min(y_px_range.end),
    );

    let area = chart.plotting_area().strip_coord_spec();
    let (px_left, py_top) = chart.as_coord_spec().translate(&(x_min, y_max));
    let (px_right, py_bottom) = chart.as_coord_spec().translate(&(x_max, y_min));
    let (ly_top, ly_bottom) = (py_top - y_offset, py_bottom - y_offset);
    let (lx_left, lx_right) = (px_left - x_offset, px_right - x_offset);

    // Top edge: translate each x tick at y_max to get its pixel position
    area.draw(&PathElement::new(
        vec![(lx_left, ly_top), (lx_right, ly_top)],
        style,
    ))?;
    for &x in x_ticks {
        let (px, py_top) = chart.as_coord_spec().translate(&(x, y_max));
        let (lx, ly) = (px - x_offset, py_top - y_offset);
        area.draw(&PathElement::new(
            vec![(lx, ly), (lx, ly + tick_size)],
            style,
        ))?;
    }

    // Right edge: translate each y tick at x_max to get its pixel position
    area.draw(&PathElement::new(
        vec![(lx_right, ly_top), (lx_right, ly_bottom)],
        style,
    ))?;
    for &y in y_ticks {
        let (px_right, py) = chart.as_coord_spec().translate(&(x_max, y));
        let (lx, ly) = (px_right - x_offset, py - y_offset);
        area.draw(&PathElement::new(
            vec![(lx, ly), (lx - tick_size, ly)],
            style,
        ))?;
    }

    // Top-right corner tick dot (optional, closes the box)
    let (px_right, py_top) = chart.as_coord_spec().translate(&(x_max, y_max));
    let (lx, ly) = (px_right - x_offset, py_top - y_offset);
    area.draw(&PathElement::new(
        vec![(lx, ly), (lx - tick_size, ly)],
        style,
    ))?;
    area.draw(&PathElement::new(
        vec![(lx, ly), (lx, ly + tick_size)],
        style,
    ))?;

    Ok(())
}

pub fn instability_plot() {
    let stable = MCMCOutput::load("data/512_x_100k_stable.mcmc").unwrap();
    let unstable = MCMCOutput::load("data/512_x_100k_unstable.mcmc").unwrap();

    let mut stable_params = stable.best_params;
    let mut unstable_params = unstable.best_params;

    (stable_params[1], stable_params[0]) = m200_c200_to_rs_rhos(stable_params[0], stable_params[1]);
    (unstable_params[1], unstable_params[0]) =
        m200_c200_to_rs_rhos(unstable_params[0], unstable_params[1]);

    instability_showcase(
        &vec![stable_params, unstable_params],
        (3e3, 1e8),
        &relhic_temperature,
        &relhic_neutral_fraction,
        ("sans-serif", 30).into_font(),
    );
}

pub fn cdm_vs_sidm_fit_plot(data: &fitting::Data) {
    let sidm_output = MCMCOutput::load("data/512_x_100k_stable.mcmc").unwrap();
    dbg!(&sidm_output.best_params);
    let (sidm_m200, sidm_c200) = (sidm_output.best_params[0], sidm_output.best_params[1]);
    let (sidm_rs, sidm_rhos) = m200_c200_to_rs_rhos(sidm_m200, sidm_c200);
    let sidm_halo = Halo::NFW(sidm_rhos, sidm_rs);

    let cdm_params = MCMCOutput::load("data/512_x_100k_sigma=0_stable.mcmc")
        .unwrap()
        .best_params;
    dbg!(&cdm_params);
    let (cdm_m200, cdm_c200) = (cdm_params[0], cdm_params[1]);
    let (cdm_rs, cdm_rhos) = m200_c200_to_rs_rhos(cdm_m200, cdm_c200);
    let cdm_halo = Halo::NFW(cdm_rhos, cdm_rs);

    let collapse_params = sidm_output
        .chain
        .into_iter()
        .zip(sidm_output.log_likelihoods)
        .filter(|(p, _)| p[2] > 0.95)
        .max_by(|(_, la), (_, lb)| la.partial_cmp(lb).unwrap())
        .map(|(p, _)| p)
        .unwrap();
    dbg!(&collapse_params);
    let (collapse_m200, collapse_c200) = (collapse_params[0], collapse_params[1]);
    let (collapse_rs, collapse_rhos) = m200_c200_to_rs_rhos(collapse_m200, collapse_c200);
    let collapse_halo = Halo::NFW(collapse_rhos, collapse_rs);

    let r_max = 1e2
        * cdm_halo
            .r_crit()
            .max(sidm_halo.r_crit())
            .max(collapse_halo.r_crit()); //shared r_max so they share r_grid

    let cdm_fit = core_collapse_background(
        relhic_temperature,
        relhic_neutral_fraction,
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
        relhic_neutral_fraction,
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
        relhic_neutral_fraction,
        collapse_rhos,
        collapse_rs,
        collapse_tau,
        Some(collapse_params[3]),
        (INNER_BOUND, r_max),
        false,
    );

    let root = SVGBackend::new("figures/cdm_vs_sidm_stable.svg", (1024, 768)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    draw_functions_on_area(
        &root,
        &cdm_fit.0,
        &vec![cdm_fit.1, sidm_fit.1, collapse_fit.1],
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
    .unwrap();

    // Build Inset
    println!("Building Inset!");

    let r_dm_min = 1e-2 * cdm_rs.min(sidm_rs).min(collapse_rs);
    let r_dm_max = 1e2 * cdm_rs.max(sidm_rs).max(collapse_rs);

    let inset = root.clone().shrink((170, 290), (512, 400));
    inset.fill(&WHITE).unwrap();

    let rho_dm_max = cdm_rhos.min(sidm_rhos).min(collapse_rhos) * 1e3;
    let rho_dm_min = cdm_rhos.max(sidm_rhos).max(collapse_rhos) * 1e-3;

    let mut inset_chart = ChartBuilder::on(&inset)
        .margin(5)
        .x_label_area_size(60)
        .y_label_area_size(120)
        .build_cartesian_2d(
            (r_dm_min..r_dm_max).log_scale(),
            (rho_dm_min..rho_dm_max).log_scale(),
        )
        .unwrap();
    // .set_secondary_coord(0..1, 0..1); // no-op placeholder if you don't need secondary axis

    inset_chart
        .configure_mesh()
        .x_desc("r (kpc)")
        .y_desc("Dark Matter Density (Mₛᵤₙ kpc⁻³)")
        .x_label_style(("sans-serif", 25).into_font())
        .y_label_style(("sans-serif", 25).into_font())
        .x_labels(3)
        .y_labels(3)
        .x_label_formatter(&fmt_num)
        .y_label_formatter(&fmt_num)
        .draw()
        .unwrap();

    // Build/Plot density profiles

    const R_GRID_NUM: usize = 1000;
    let r_range: Vec<f64> = (0..R_GRID_NUM)
        .into_iter()
        .map(|i: usize| -> f64 {
            (r_dm_min.ln() + (i as f64 / (R_GRID_NUM - 1) as f64) * (r_dm_max.ln() - r_dm_min.ln()))
                .exp()
        })
        .collect();

    for (rs, rhos, m200, c200, tau, color) in [
        (cdm_rs, cdm_rhos, cdm_m200, cdm_c200, 0.0, &BLACK),
        (
            sidm_rs,
            sidm_rhos,
            sidm_m200,
            sidm_c200,
            sidm_output.best_params[2],
            &BLUE,
        ),
        (
            collapse_rs,
            collapse_rhos,
            collapse_m200,
            collapse_c200,
            collapse_params[2],
            &RED,
        ),
    ] {
        let func = parametic_core_collapse(rs, rhos, tau);
        let dm_rho_pts: Vec<(f64, f64)> = r_range.iter().map(|&r| (r, func(r))).collect();
        inset_chart
            .draw_series(LineSeries::new(dm_rho_pts, color))
            .unwrap()
            .label(format!("M₂₀₀ = {:.2e}, C₂₀₀= {:.2}", m200, c200))
            .legend(move |(x, y)| {
                PathElement::new(
                    vec![(x, y), (x + 20, y)],
                    ShapeStyle {
                        color: color.mix(1.0),
                        filled: false,
                        stroke_width: 2,
                    },
                )
            });
    }

    // Configure and draw legend
    inset_chart
        .configure_series_labels()
        .label_font(("sans-serif", 25).into_font())
        .position(SeriesLabelPosition::UpperRight)
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    let (x_ticks, y_ticks) = (
        log_ticks(r_dm_min, r_dm_max, 0),
        log_ticks(rho_dm_min, rho_dm_max, 0),
    );

    draw_ticks_top_and_right(
        &inset_chart,
        &x_ticks,
        &y_ticks,
        (r_dm_min, r_dm_max),
        (rho_dm_min, rho_dm_max),
        5,
        BLACK.stroke_width(1),
    )
    .unwrap();

    root.present().unwrap();

    svg_to_pdf("figures/cdm_vs_sidm_stable").unwrap();
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

        if let Ok(output) = MCMCOutput::load(&(String::from("data/") + &file_name + "_stable.mcmc")) {
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
        &(String::from("figures/cross_section_vs_deviation_stable.svg")),
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
