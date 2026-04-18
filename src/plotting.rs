use ensemble_mcmc::MCMCOutput;
use plotters::prelude::*;

use crate::constants::*;
use crate::corner_plot::*;
use crate::halo::{McrSource, deviation, m200_c200_to_rs_rhos};
// use crate::logging::load_file;
use crate::utils::make_file_name;
use std::f64::consts::PI;

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
    //let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
    let root = SVGBackend::new(filename, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;

    let data_spread = match data {
        Some(data_points) => Some(data_points.last().unwrap().0 - data_points[0].0),
        None => None,
    };

    let x_range = match data {
        Some(data_points) => 0.5 * data_points[0].0..2.0 * data_points.last().unwrap().0,
        // Some(data_points) => {
        //     (data_points[0].0 - 0.1 * data_spread.unwrap())
        //         ..(data_points.last().unwrap().0 + 0.1 * data_spread.unwrap())
        // }
        None => x_points[0]..x_points[x_points.len() - 1],
    };

    // let x_range = x_points[0]..x_points[x_points.len() - 1];

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

    //println!("y_min = {:.3}, y_max={:.3}", y_min, y_max);

    let y_range = (y_min + 1e-4)
        * match y_min.signum() {
            1.0 => 0.9,
            -1.0 => 1.1,
            _ => panic!("number has no sign, is probably NaN"),
        }
        ..y_max
            * match y_max.signum() {
                1.0 => 3.0,
                -1.0 => 0.9,
                _ => panic!("number has no sign, is probably NaN"),
            };

    let log_scale = false;

    let x_range = x_range.log_scale();
    let y_range = y_range.log_scale();
    let log_scale = true;

    let mut chart = ChartBuilder::on(&root)
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

    root.present()?;
    //println!("Plot saved as {}", filename);
    Ok(())
}

pub fn create_cross_section_deviation_relation_plot(
    num_walkers: usize,
    steps: usize,
    prior: &crate::fitting::Prior,
    cross_sections: &[f64],
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut data: Vec<(f64, f64)> = Vec::with_capacity(cross_sections.len());
    let mut y_err: Vec<(f64, f64)> = Vec::with_capacity(cross_sections.len());

    for cross_section in cross_sections {
        let file_name = make_file_name(num_walkers, steps, prior, &Some(*cross_section));

        if let Ok(output) = MCMCOutput::load(&(String::from("data/") + &file_name + ".mcmc")) {
            let chain = output.chain;

            let mut deviation_chain = Vec::with_capacity(chain.len());
            let mut mean_deviation = 0.0;
            let mut mean_square_deviation = 0.0;
            for params in &chain {
                let deviation = deviation(params[0], params[1], McrSource::DiemerJoyce2019);

                mean_deviation += deviation;
                mean_square_deviation += deviation.powi(2);
                deviation_chain.push(deviation);
            }

            mean_deviation /= chain.len() as f64;
            mean_square_deviation /= chain.len() as f64;
            let var_deviation = mean_square_deviation - mean_deviation.powi(2);

            data.push((*cross_section, mean_deviation));
            y_err.push((
                mean_deviation - var_deviation.sqrt(),
                mean_deviation + var_deviation.sqrt(),
            ));
        } else {
            println!(
                "Unable to find data for cross_section = {cross_section}! Has this MCMC been run?"
            )
        }
    }

    plot_functions(
        &Vec::new(),
        &Vec::new(),
        &(String::from("figures/cross_section_vs_deviation.svg")),
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

    plot_histogram(
        &root.margin(5, 5 + X_LABEL_HEIGHT, 5 + Y_LABEL_WIDTH, 5),
        &deviation_chain,
        "deviation from cosmological median",
        false,
        (0, 0),
        1,
        &bounds,
        &marked_values,
        None,
        font,
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

    plot_histogram(
        &root.margin(5, 5 + X_LABEL_HEIGHT, 5 + Y_LABEL_WIDTH, 5),
        &t_sigma_m_chain,
        "log(t sigma/m)",
        false,
        (0, 0),
        1,
        &bounds,
        &marked_values,
        None, //Some(&[&prob_dens]),
        font,
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
