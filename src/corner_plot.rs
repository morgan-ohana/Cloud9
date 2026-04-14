use plotters::coord::Shift;
use plotters::coord::ranged1d::ValueFormatter;
use plotters::element::DashedPathElement;
use plotters::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

//Paper
// const LABEL_WIDTH: u32 = 30;
// const Y_LABEL_PAD: u32 = 15;
// const GAP: u32 = 5;
//Presentation
pub const X_LABEL_HEIGHT: u32 = 60;
pub const Y_LABEL_WIDTH: u32 = 100;
pub const MARGIN: u32 = 10;
pub const GAP: u32 = 15;

pub fn fmt_num(num: &f64) -> String {
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

pub fn create_corner_plot(
    chain: &[Vec<f64>],
    marked_points: &[&[f64]],
    param_names: &[&str],
    log_scales: &[bool],
    output_path: &str,
    bounds: &[[f64; 2]],
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Extract parameter columns
    let num_params = chain[0].len();

    let mut params: Vec<Vec<f64>> = vec![Vec::new(); num_params];
    for point in chain {
        for i in 0..num_params {
            params[i].push(point[i]);
        }
    }

    let output_path = String::from(output_path);
    let svg_path = &(output_path.clone() + ".svg");
    let pdf_path = &(output_path.clone() + ".pdf");
    // Create plot area

    // Plot Width = Plot Height. Demand square
    let plot_width = 1600;
    let cell_size = plot_width / num_params as u32;

    let inner_width = plot_width + Y_LABEL_WIDTH + (num_params as u32 - 1) * GAP;
    let inner_height = plot_width + X_LABEL_HEIGHT + (num_params as u32 - 1) * GAP;

    let root_width = inner_width + 2 * MARGIN;
    let root_height = inner_height + 2 * MARGIN;

    let root = SVGBackend::new(svg_path, (root_width, root_height)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_break_points: Vec<u32> = (1..num_params)
        .map(|i| Y_LABEL_WIDTH + (i as u32 * cell_size) + ((i as u32 - 1) * GAP))
        .collect();

    let y_break_points: Vec<u32> = (1..num_params)
        .map(|i| (i as u32 * cell_size) + ((i as u32 - 1) * GAP))
        .collect();

    // Split into subplots
    let sub_areas = root
        .margin(MARGIN, MARGIN, MARGIN, MARGIN)
        .split_by_breakpoints(x_break_points, y_break_points);

    // Plot each cell
    for row in 0..num_params {
        for col in 0..num_params {
            let idx = row * num_params + col;
            let drawing_area = &sub_areas[idx];

            if row == col {
                // Diagonal: Histogram
                let marked_values: Vec<f64> = marked_points.iter().map(|p| p[row]).collect();
                plot_histogram(
                    drawing_area,
                    &params[row],
                    param_names[row],
                    log_scales[row],
                    (row, col),
                    num_params,
                    &bounds[row],
                    &marked_values,
                    None,
                    font.clone(),
                )?;
            } else if row > col {
                // Lower triangle: 2D scatter/density
                let marked_2d: Vec<(f64, f64)> =
                    marked_points.iter().map(|p| (p[col], p[row])).collect();
                plot_2d_scatter(
                    drawing_area,
                    &params[col],
                    &params[row],
                    param_names,
                    log_scales,
                    (row, col),
                    num_params,
                    &bounds,
                    &marked_2d,
                    font.clone(),
                )?;
            } else {
                // Upper triangle: Correlation/contour or leave empty
                plot_correlation(drawing_area, &params[col], &params[row], font.clone())?;
            }
        }
    }

    root.present()?;
    println!("Corner plot saved to: {}.svg", output_path);

    let status = std::process::Command::new("inkscape")
        .args(&["--export-type=pdf", "--export-filename", pdf_path, svg_path])
        .status()?;

    if status.success() {
        Ok(())
    } else {
        Err(format!("Inkscape failed with status: {}", status).into())
    }
}

pub fn plot_histogram(
    area: &DrawingArea<SVGBackend, Shift>,
    data: &[f64],
    param_name: &str,
    log_scale: bool,
    (row, col): (usize, usize),
    num_params: usize,
    bounds: &[f64; 2],
    marked_values: &[f64],
    extra_functions: Option<&[&dyn Fn(f64) -> f64]>,
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Calculate bins
    let min = bounds[0];
    let max = bounds[1];

    let n_bins = 500;
    let edges = {
        match log_scale {
            true => {
                assert!(min > 0.0, "Cannot use log scale for non-positive values");
                let mut edges = vec![0.0; n_bins + 1];
                for i in 0..n_bins + 1 {
                    let t = (i as f64) / (n_bins as f64);
                    edges[i] = (min.ln() + t * (max.ln() - min.ln())).exp();
                }
                edges
            }
            false => {
                let mut edges = vec![0.0; n_bins + 1];
                for i in 0..n_bins + 1 {
                    let t = (i as f64) / (n_bins as f64);
                    edges[i] = min + t * (max - min);
                }
                edges
            }
        }
    };

    let spacing = match log_scale {
        true => (max.ln() - min.ln()) / (n_bins as f64),
        false => (max - min) / (n_bins as f64),
    };

    let mut bins: Vec<u32> = vec![0; n_bins];

    let mut counts = 0;
    for &value in data {
        if value >= max || value < min {
            continue; // skip out of range values
        }

        let bin_idx = match log_scale {
            true => ((value.ln() - min.ln()) / spacing).floor() as usize,
            false => ((value - min) / spacing).floor() as usize,
        };

        bins[bin_idx] += 1;
        counts += 1;
    }

    if counts == 0 {
        return Err(format!(
            "No data points fall within the provided bounds in {param_name} histogram"
        )
        .into());
    }

    let mut density = vec![0.0; n_bins];
    let mut max_density: f64 = 0.0;
    for i in 0..bins.len() {
        let width = edges[i + 1] - edges[i];
        density[i] = bins[i] as f64 / (counts as f64 * width);
        max_density = max_density.max(density[i]);
    }

    // Create chart for histogram
    let mut chart_builder = ChartBuilder::on(area);

    if col == 0 {
        chart_builder.y_label_area_size(Y_LABEL_WIDTH);
    } else {
        chart_builder.margin_left(GAP);
    }

    if row == num_params - 1 {
        chart_builder.x_label_area_size(X_LABEL_HEIGHT);
    } else {
        chart_builder.margin_bottom(GAP);
    }

    match log_scale {
        true => {
            let mut chart = chart_builder
                //.caption(param_name, ("sans-serif", 15).into_font())
                .build_cartesian_2d((min..max).log_scale(), 0.0..max_density * 1.1)?;

            draw_hist_content(
                &mut chart,
                param_name,
                &density,
                &max_density,
                &edges,
                font.clone(),
            )?;

            // Draw marked values with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &marked_value) in marked_values.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(PathElement::new(
                    vec![(marked_value, 0.0), (marked_value, max_density * 1.1)],
                    color,
                )))?;
            }

            // Add KDE curve
            plot_kde(&mut chart, data, min, max, log_scale)?;

            // Add extra functions
            if let Some(extra_funcs) = extra_functions {
                let x_points: Vec<f64> = (0..=1000)
                    .map(|i| (min.ln() + (i as f64) * (max.ln() - min.ln()) / 1000.0).exp())
                    .collect();
                for extra_func in extra_funcs {
                    let points: Vec<(f64, f64)> =
                        x_points.iter().map(|x| (*x, extra_func(*x))).collect();
                    let _series = chart.draw_series(LineSeries::new(points, &GREEN));
                }
            }

            // Draw Legend
            chart
                .configure_series_labels()
                .label_font(font)
                .position(SeriesLabelPosition::UpperRight)
                .background_style(&WHITE.mix(0.8))
                .border_style(&BLACK)
                .draw()?;
        }
        false => {
            let mut chart = chart_builder
                //.caption(param_name, ("sans-serif", 15).into_font())
                .build_cartesian_2d(min..max, 0.0..max_density * 1.1)?;

            draw_hist_content(
                &mut chart,
                param_name,
                &density,
                &max_density,
                &edges,
                font.clone(),
            )?;

            // Draw marked values with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &marked_value) in marked_values.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(PathElement::new(
                    vec![(marked_value, 0.0), (marked_value, max_density * 1.1)],
                    color,
                )))?;
            }

            // Add KDE curve
            plot_kde(&mut chart, data, min, max, log_scale)?;

            // Add extra functions
            if let Some(extra_funcs) = extra_functions {
                let x_points: Vec<f64> = (0..=1000)
                    .map(|i| min + (i as f64) * (max - min) / 1000.0)
                    .collect();
                for extra_func in extra_funcs {
                    let points: Vec<(f64, f64)> =
                        x_points.iter().map(|x| (*x, extra_func(*x))).collect();
                    let _series = chart.draw_series(LineSeries::new(points, &GREEN));
                }
            }

            // Draw Legend
            chart
                .configure_series_labels()
                .label_font(font)
                .position(SeriesLabelPosition::UpperRight)
                .background_style(&WHITE.mix(0.8))
                .border_style(&BLACK)
                .draw()?;
        }
    };

    Ok(())
}

fn draw_hist_content<
    X: Ranged<ValueType = f64> + ValueFormatter<f64>,
    Y: Ranged<ValueType = f64> + ValueFormatter<f64>,
>(
    chart: &mut ChartContext<SVGBackend, Cartesian2d<X, Y>>,
    param_name: &str,
    density: &Vec<f64>,
    max_density: &f64,
    edges: &Vec<f64>,
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    chart
        .configure_mesh()
        .x_label_style(font.clone())
        .x_desc(param_name) // X-axis label
        .y_label_style(font.clone())
        .y_desc("Probability Density") // Y-axis label
        .x_label_formatter(&fmt_num)
        .y_label_formatter(&fmt_num)
        .label_style(font.clone())
        .draw()?;

    // Plot histogram bars
    let mut bins = Vec::with_capacity(density.len());
    for i in 0..density.len() {
        let count = density[i];

        chart.draw_series(std::iter::once(Rectangle::new(
            [(edges[i], 0.0), (edges[i + 1], count)],
            BLUE.mix(0.5).filled(),
        )))?;

        bins.push((i, density[i]));
    }

    bins.sort_by(|bin_a, bin_b| bin_b.1.partial_cmp(&bin_a.1).unwrap());

    let mode_idx = bins[0].0;
    let mode_point = 0.5 * (edges[mode_idx] + edges[mode_idx + 1]);
    println!("Mode {param_name} = {:.5e}", mode_point);

    let mut cumulative_probability = 0.0;
    let sigma_levels = [0.6827, 0.9545, 0.9973];
    let mut included_indices = vec![Vec::new(); sigma_levels.len()];
    for i in 0..bins.len() {
        let idx = bins[i].0;
        cumulative_probability += density[idx] * (edges[idx + 1] - edges[idx]);
        for n in 0..sigma_levels.len() {
            if cumulative_probability > sigma_levels[n] {
                continue;
            }

            included_indices[n].push(idx);
        }

        if &cumulative_probability > sigma_levels.last().unwrap() {
            break;
        }
    }

    let mut interval_edges = Vec::new();
    for n in 0..sigma_levels.len() {
        included_indices[n].sort();

        let (left_idx, right_idx): (usize, usize) =
            (included_indices[n][0], *included_indices[n].last().unwrap());

        let mut included_prob = 0.0;
        for i in left_idx..=right_idx {
            included_prob += density[i] * (edges[i + 1] - edges[i]);
        }
        println!("{} = {:.4}", sigma_levels[n], included_prob);

        // get positions of edges that bracket 1-sigma interval
        let left_edge = 0.5 * (edges[left_idx] + edges[left_idx + 1]);
        let right_edge = 0.5 * (edges[right_idx] + edges[right_idx + 1]);

        println!(
            "{} sigma bounds: ({:.5}, {:.5})",
            n + 1,
            left_edge,
            right_edge
        );

        interval_edges.push((left_edge, right_edge));
    }

    chart
        .draw_series(std::iter::once(PathElement::new(
            vec![(mode_point, 0.0), (mode_point, max_density * 1.1)],
            &BLACK,
        )))?
        .label(format!(
            "{param_name} = {}  ({}, {})",
            fmt_num(&mode_point),
            fmt_num(&interval_edges[0].0),
            fmt_num(&interval_edges[0].1)
        ))
        .legend(move |(x, y)| {
            PathElement::new(
                vec![(x, y), (x + 20, y)],
                ShapeStyle {
                    color: BLACK.mix(1.0),
                    filled: false,
                    stroke_width: 2,
                },
            )
        });

    chart.draw_series(std::iter::once(DashedPathElement::new(
        vec![
            (interval_edges[0].0, 0.0),
            (interval_edges[0].0, max_density * 1.1),
        ],
        10,
        10,
        &BLACK,
    )))?;

    chart.draw_series(std::iter::once(DashedPathElement::new(
        vec![
            (interval_edges[0].1, 0.0),
            (interval_edges[0].1, max_density * 1.1),
        ],
        10,
        10,
        &BLACK,
    )))?;

    Ok(())
}

fn plot_2d_scatter(
    area: &DrawingArea<SVGBackend, Shift>,
    x_data: &[f64],
    y_data: &[f64],
    param_names: &[&str],
    log_scales: &[bool],
    (row, col): (usize, usize),
    num_params: usize,
    bounds: &[[f64; 2]],
    marked_points_2d: &[(f64, f64)],
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    let [x_min, x_max] = bounds[col];
    let [y_min, y_max] = bounds[row];

    let mut chart_builder = ChartBuilder::on(area);

    if col == 0 {
        chart_builder.y_label_area_size(Y_LABEL_WIDTH);
    } else {
        chart_builder.margin_left(GAP);
    }

    if row == num_params - 1 {
        chart_builder.x_label_area_size(X_LABEL_HEIGHT);
    } else {
        chart_builder.margin_bottom(GAP);
    }

    match (log_scales[col], log_scales[row]) {
        (true, true) => {
            let mut chart = chart_builder
                .build_cartesian_2d((x_min..x_max).log_scale(), (y_min..y_max).log_scale())?;

            draw_scatter_content(&mut chart, param_names, (row, col), font)?;

            // Draw marked points with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &(marked_x, marked_y)) in marked_points_2d.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(Circle::new(
                    (marked_x, marked_y),
                    5,
                    color.filled(),
                )))?;
            }
        }
        (true, false) => {
            let mut chart =
                chart_builder.build_cartesian_2d((x_min..x_max).log_scale(), y_min..y_max)?;

            draw_scatter_content(&mut chart, param_names, (row, col), font)?;

            // Draw marked points with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &(marked_x, marked_y)) in marked_points_2d.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(Circle::new(
                    (marked_x, marked_y),
                    5,
                    color.filled(),
                )))?;
            }
        }
        (false, true) => {
            let mut chart =
                chart_builder.build_cartesian_2d(x_min..x_max, (y_min..y_max).log_scale())?;

            draw_scatter_content(&mut chart, param_names, (row, col), font)?;

            // Draw marked points with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &(marked_x, marked_y)) in marked_points_2d.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(Circle::new(
                    (marked_x, marked_y),
                    5,
                    color.filled(),
                )))?;
            }
        }
        (false, false) => {
            let mut chart = chart_builder.build_cartesian_2d(x_min..x_max, y_min..y_max)?;

            draw_scatter_content(&mut chart, param_names, (row, col), font)?;

            // Draw marked points with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &(marked_x, marked_y)) in marked_points_2d.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(Circle::new(
                    (marked_x, marked_y),
                    5,
                    color.filled(),
                )))?;
            }
        }
    };

    // Add contour lines for density
    plot_2d_contours(
        area,
        x_data,
        y_data,
        x_min,
        x_max,
        y_min,
        y_max,
        log_scales,
        param_names,
        (row, col),
        num_params,
    )?;

    Ok(())
}

fn draw_scatter_content<
    X: Ranged<ValueType = f64> + ValueFormatter<f64>,
    Y: Ranged<ValueType = f64> + ValueFormatter<f64>,
>(
    chart: &mut ChartContext<SVGBackend, Cartesian2d<X, Y>>,
    param_names: &[&str],
    (row, col): (usize, usize),
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut mesh = chart.configure_mesh();

    mesh.x_label_style(font.clone())
        .x_desc(param_names[col]) // X-axis label
        .y_label_style(font.clone())
        .y_desc(param_names[row]) // Y-axis label
        .x_label_formatter(&fmt_num)
        .y_label_formatter(&fmt_num)
        .label_style(font.clone());

    mesh.draw()?;

    Ok(())
}

fn plot_kde<
    X: Ranged<ValueType = f64> + ValueFormatter<f64>,
    Y: Ranged<ValueType = f64> + ValueFormatter<f64>,
>(
    chart: &mut ChartContext<SVGBackend, Cartesian2d<X, Y>>,
    data: &[f64],
    min: f64,
    max: f64,
    log_scale: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Simple KDE using Gaussian kernel

    let n = data.len() as f64;

    let bandwidth = if log_scale {
        let log_data: Vec<f64> = data.iter().map(|x| x.ln()).collect();
        let mean = log_data.iter().sum::<f64>() / n;
        let var = log_data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        let sigma = var.sqrt();
        1.06 * sigma * n.powf(-0.2)
    } else {
        let mean = data.iter().sum::<f64>() / n;
        let var = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
        let sigma = var.sqrt();
        1.06 * sigma * n.powf(-0.2)
    };

    // 1. Precalculate logarithms once to avoid redundant computations
    let log_data: Option<Vec<f64>> = if log_scale {
        Some(data.iter().map(|x| x.ln()).collect())
    } else {
        None
    };

    let n_points = 200;

    // Pre-calculate the normalization factor
    let norm_factor = data.len() as f64 * bandwidth * (2.0 * std::f64::consts::PI).sqrt();

    let kde_points: Vec<(f64, f64)> = (0..n_points)
        .into_par_iter()
        .map(|i| {
            let x = match log_scale {
                true => (min.ln() + (max.ln() - min.ln()) * i as f64 / (n_points - 1) as f64).exp(),
                false => min + (max - min) * i as f64 / (n_points - 1) as f64,
            };

            let mut density = 0.0;

            if log_scale {
                // 2. Compute x.ln() exactly once per KDE point
                let ln_x = x.ln();
                for &ln_point in log_data.as_ref().unwrap() {
                    let diff = (ln_x - ln_point) / bandwidth;
                    density += (-0.5 * diff * diff).exp();
                }
                density /= norm_factor * x; // Adjust for log scale
            } else {
                for &point in data {
                    let diff = (x - point) / bandwidth;
                    density += (-0.5 * diff * diff).exp();
                }
                density /= norm_factor;
            }

            (x, density)
        })
        .collect();

    chart.draw_series(LineSeries::new(kde_points, &RED))?;

    Ok(())
}

fn gaussian_kernel(sigma: f64) -> Vec<f64> {
    let radius = (3.0 * sigma).ceil() as isize;
    let mut kernel = Vec::new();
    let mut sum = 0.0;

    for i in -radius..=radius {
        let x = i as f64;
        let v = (-x * x / (2.0 * sigma * sigma)).exp();
        kernel.push(v);
        sum += v;
    }

    for v in kernel.iter_mut() {
        *v /= sum;
    }

    kernel
}

fn blur_x(field: &Vec<Vec<f64>>, sigma: f64) -> Vec<Vec<f64>> {
    let kernel = gaussian_kernel(sigma);
    let r = (kernel.len() / 2) as isize;
    let n = field.len();
    let m = field[0].len();

    let mut out = vec![vec![0.0; m]; n];

    for i in 0..n {
        for j in 0..m {
            let mut sum = 0.0;
            for k in -r..=r {
                let jj = (j as isize + k).clamp(0, (m - 1) as isize) as usize;
                sum += field[i][jj] * kernel[(k + r) as usize];
            }
            out[i][j] = sum;
        }
    }
    out
}

fn blur_y(field: &Vec<Vec<f64>>, sigma: f64) -> Vec<Vec<f64>> {
    let kernel = gaussian_kernel(sigma);
    let r = (kernel.len() / 2) as isize;
    let n = field.len();
    let m = field[0].len();

    let mut out = vec![vec![0.0; m]; n];

    for i in 0..n {
        for j in 0..m {
            let mut sum = 0.0;
            for k in -r..=r {
                let ii = (i as isize + k).clamp(0, (n - 1) as isize) as usize;
                sum += field[ii][j] * kernel[(k + r) as usize];
            }
            out[i][j] = sum;
        }
    }
    out
}

fn gaussian_smooth(field: &Vec<Vec<f64>>, sigma: f64) -> Vec<Vec<f64>> {
    let tmp = blur_x(field, sigma);
    blur_y(&tmp, sigma)
}

fn plot_2d_contours(
    area: &DrawingArea<SVGBackend, Shift>,
    x_data: &[f64],
    y_data: &[f64],
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    log_scales: &[bool],
    param_names: &[&str],
    (row, col): (usize, usize),
    num_params: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    // Simple 2D density estimation
    let grid_size = 100;

    let mut edges = vec![vec![(0.0, 0.0); grid_size + 1]; grid_size + 1];
    for i in 0..=grid_size {
        for j in 0..=grid_size {
            let x_edge = match log_scales[col] {
                true => {
                    (x_min.ln() + (x_max.ln() - x_min.ln()) * i as f64 / grid_size as f64).exp()
                }
                false => x_min + (x_max - x_min) * i as f64 / grid_size as f64,
            };
            let y_edge = match log_scales[row] {
                true => {
                    (y_min.ln() + (y_max.ln() - y_min.ln()) * j as f64 / grid_size as f64).exp()
                }
                false => y_min + (y_max - y_min) * j as f64 / grid_size as f64,
            };
            edges[i][j] = (x_edge, y_edge);
        }
    }

    let mut density = vec![vec![0.0; grid_size]; grid_size];
    let mut count = 0;
    for (&x, &y) in x_data.iter().zip(y_data.iter()) {
        if x >= x_max || x < x_min || y >= y_max || y < y_min {
            continue;
        }

        let x_bin = match log_scales[col] {
            true => ((x.ln() - x_min.ln()) / (x_max.ln() - x_min.ln()) * grid_size as f64).floor()
                as usize,
            false => ((x - x_min) / (x_max - x_min) * grid_size as f64).floor() as usize,
        };
        let y_bin = match log_scales[row] {
            true => ((y.ln() - y_min.ln()) / (y_max.ln() - y_min.ln()) * grid_size as f64).floor()
                as usize,
            false => ((y - y_min) / (y_max - y_min) * grid_size as f64).floor() as usize,
        };

        count += 1;
        density[x_bin][y_bin] += 1.0;
    }

    if count == 0 {
        return Err(format!(
            "No data points fall within the provided bounds in {} x {} contour plot",
            param_names[row], param_names[col]
        )
        .into());
    }

    //Normalize
    for i in 0..grid_size {
        for j in 0..grid_size {
            let area = (edges[i + 1][j].0 - edges[i][j].0) * (edges[i][j + 1].1 - edges[i][j].1);
            density[i][j] /= area * count as f64
        }
    }

    let y_margin_width = match col {
        0 => Y_LABEL_WIDTH,
        _ => GAP,
    };
    let x_margin_width = if row == num_params - 1 {
        X_LABEL_HEIGHT
    } else {
        GAP
    };

    match (log_scales[col], log_scales[row]) {
        (true, true) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(x_margin_width)
                .build_cartesian_2d((x_min..x_max).log_scale(), (y_min..y_max).log_scale())?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
        (true, false) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(x_margin_width)
                .build_cartesian_2d((x_min..x_max).log_scale(), y_min..y_max)?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
        (false, true) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(x_margin_width)
                .build_cartesian_2d(x_min..x_max, (y_min..y_max).log_scale())?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
        (false, false) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(x_margin_width)
                .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
    }

    Ok(())
}

fn draw_contour_content(
    chart: &mut ChartContext<
        SVGBackend,
        Cartesian2d<impl Ranged<ValueType = f64>, impl Ranged<ValueType = f64>>,
    >,
    density: &Vec<Vec<f64>>,
    edges: &Vec<Vec<(f64, f64)>>,
    grid_size: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let density = gaussian_smooth(&density, 1.0);

    let mut cells = Vec::with_capacity(grid_size.pow(2));
    for i in 0..grid_size {
        for j in 0..grid_size {
            let x_start = edges[i][j].0;
            let x_end = edges[i + 1][j].0;
            let y_start = edges[i][j].1;
            let y_end = edges[i][j + 1].1;

            cells.push((density[i][j], (x_start, y_start), (x_end, y_end), (i, j)));
        }
    }

    cells.sort_by(|cell_a, cell_b| cell_b.0.partial_cmp(&cell_a.0).unwrap());

    const CONTOUR_LEVELS: [f64; 5] = [0.118, 0.393, 0.675, 0.864, 0.956]; //[0.5, 1, 1.5, 2, 2.5] sigma for 2D gaussian
    let contour_colors = [
        RGBColor(255, 100, 100),
        RGBColor(255, 200, 100),
        RGBColor(100, 255, 200),
        RGBColor(100, 200, 255),
        RGBColor(100, 100, 255),
    ];
    let mut cumulative_probabilities = vec![vec![0.0; grid_size]; grid_size];
    let mut cumulative_probability = 0.0;
    let mut levels = vec![vec![CONTOUR_LEVELS.len(); grid_size]; grid_size];
    let mut density_levels = [1.0; 5];

    for cell in cells {
        let area = (edges[cell.3.0 + 1][cell.3.1].0 - edges[cell.3.0][cell.3.1].0)
            * (edges[cell.3.0][cell.3.1 + 1].1 - edges[cell.3.0][cell.3.1].1);
        cumulative_probability += cell.0 * area;
        cumulative_probabilities[cell.3.0][cell.3.1] = cumulative_probability;

        let mut color_opt = None;
        for i in 0..CONTOUR_LEVELS.len() {
            if cumulative_probability <= CONTOUR_LEVELS[i] {
                levels[cell.3.0][cell.3.1] = i;
                density_levels[i] = cell.0;
                color_opt = Some(contour_colors[i]);
                break;
            }
        }

        if let Some(color) = color_opt {
            chart.draw_series(std::iter::once(Rectangle::new(
                [cell.1, cell.2],
                color.mix(0.3).filled(),
            )))?;
        }
    }

    // Draw contours
    draw_contours(chart, density_levels[1], &density, edges, grid_size)?;
    draw_contours(chart, density_levels[3], &density, edges, grid_size)?;

    Ok(())
}

#[derive(Clone, Copy)]
struct Pt(f64, f64);

impl PartialEq for Pt {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).hypot(self.1 - other.1) < 1e-9
    }
}

impl Pt {
    fn normalize(&self, x_bounds: (f64, f64), y_bounds: (f64, f64)) -> Self {
        Pt(
            (self.0 - x_bounds.0) / (x_bounds.1 - x_bounds.0),
            (self.1 - y_bounds.0) / (y_bounds.1 - y_bounds.0),
        )
    }
}

fn interpolate(p1: Pt, p2: Pt, v1: f64, v2: f64, level: f64) -> Pt {
    let t = (level - v1) / (v2 - v1);
    Pt(p1.0 + t * (p2.0 - p1.0), p1.1 + t * (p2.1 - p1.1))
}

fn cell_to_vertex(levels: &Vec<Vec<f64>>, grid_size: usize) -> Vec<Vec<f64>> {
    let mut v = vec![vec![0.0; grid_size + 1]; grid_size + 1];

    // Intentionally avoid looping over edge verticies to leave them at 0
    for i in 0..=grid_size {
        for j in 0..=grid_size {
            let mut sum = 0.0;
            let mut count: f64 = 0.0;

            // ---------------------
            // |         |         |
            // | (i-1,j) |  (i,j)  |
            // |         |         |
            // --------(i,j)--------
            // |         |         |
            // |(i-1,j-1)| (i,j-1) |
            // |         |         |
            // ---------------------

            for di in [-1, 0] {
                for dj in [-1, 0] {
                    let ci = i as isize + di;
                    let cj = j as isize + dj;
                    if ci >= 0 && cj >= 0 && ci < grid_size as isize && cj < grid_size as isize {
                        sum += levels[ci as usize][cj as usize];
                        count += 1.0;
                    }
                }
            }

            v[i][j] = sum / count;
        }
    }

    v
}

fn marching_squares(
    levels: &Vec<Vec<f64>>,
    edges: &Vec<Vec<(f64, f64)>>,
    level: f64,
    grid_size: usize,
) -> Vec<(Pt, Pt)> {
    let mut segments = Vec::new();
    let level = level as f64;

    for i in 0..grid_size {
        for j in 0..grid_size {
            let v00 = levels[i][j];
            let v10 = levels[i + 1][j];
            let v01 = levels[i][j + 1];
            let v11 = levels[i + 1][j + 1];

            let mut mask = 0;
            if v00 > level {
                mask |= 1;
            }
            if v10 > level {
                mask |= 2;
            }
            if v11 > level {
                mask |= 4;
            }
            if v01 > level {
                mask |= 8;
            }

            if mask == 0 || mask == 15 {
                continue;
            }

            let p00 = Pt(edges[i][j].0, edges[i][j].1);
            let p10 = Pt(edges[i + 1][j].0, edges[i + 1][j].1);
            let p01 = Pt(edges[i][j + 1].0, edges[i][j + 1].1);
            let p11 = Pt(edges[i + 1][j + 1].0, edges[i + 1][j + 1].1);

            let mut edge_points = Vec::new();

            if (v00 > level) != (v10 > level) {
                edge_points.push(interpolate(p00, p10, v00, v10, level));
            }
            if (v10 > level) != (v11 > level) {
                edge_points.push(interpolate(p10, p11, v10, v11, level));
            }
            if (v11 > level) != (v01 > level) {
                edge_points.push(interpolate(p11, p01, v11, v01, level));
            }
            if (v01 > level) != (v00 > level) {
                edge_points.push(interpolate(p01, p00, v01, v00, level));
            }

            if mask == 5 || mask == 10 {
                // Use asymptotic decider
                let center = (v00 + v10 + v01 + v11) / 4.0;
                if center > level {
                    // connect high-valued corners
                    if mask == 5 {
                        edge_points = vec![
                            interpolate(p00, p01, v00, v01, level),
                            interpolate(p10, p11, v10, v11, level),
                        ];
                    } else {
                        edge_points = vec![
                            interpolate(p00, p10, v00, v10, level),
                            interpolate(p01, p11, v01, v11, level),
                        ];
                    }
                } else {
                    // connect low-valued corners
                    if mask == 5 {
                        edge_points = vec![
                            interpolate(p00, p10, v00, v10, level),
                            interpolate(p01, p11, v01, v11, level),
                        ];
                    } else {
                        edge_points = vec![
                            interpolate(p00, p01, v00, v01, level),
                            interpolate(p10, p11, v10, v11, level),
                        ];
                    }
                }
            }

            segments.push((edge_points[0], edge_points[1]));
        }
    }

    segments
}

fn stitch_segments(
    segments: &Vec<(Pt, Pt)>,
    x_bounds: (f64, f64),
    y_bounds: (f64, f64),
) -> Vec<Vec<Pt>> {
    let mut loops = Vec::new();
    let mut unused = segments.clone();

    while let Some((start, end)) = unused.pop() {
        let mut loop_pts = vec![start, end];

        loop {
            let last = *loop_pts.last().unwrap();
            let norm_last = last.normalize(x_bounds, y_bounds);
            let mut found = None;

            for (k, (a, b)) in unused.iter().enumerate() {
                let norm_a = a.normalize(x_bounds, y_bounds);
                let norm_b = b.normalize(x_bounds, y_bounds);

                if norm_a == norm_last {
                    found = Some((k, *b));
                    break;
                } else if norm_b == norm_last {
                    found = Some((k, *a));
                    break;
                }
            }

            if let Some((idx, next)) = found {
                unused.swap_remove(idx);
                if next == loop_pts[0] {
                    // Close the loop
                    loop_pts.push(next);
                    break;
                }
                loop_pts.push(next);
            } else {
                break;
            }
        }

        // Extend backward from the start (handles the case where the seed segment
        // was picked from the middle of an open chain)
        loop {
            let first = loop_pts[0];
            let norm_first = first.normalize(x_bounds, y_bounds);
            let mut found = None;

            for (k, (a, b)) in unused.iter().enumerate() {
                let norm_a = a.normalize(x_bounds, y_bounds);
                let norm_b = b.normalize(x_bounds, y_bounds);

                if norm_a == norm_first {
                    found = Some((k, *b));
                    break;
                } else if norm_b == norm_first {
                    found = Some((k, *a));
                    break;
                }
            }

            if let Some((idx, prev)) = found {
                unused.swap_remove(idx);
                loop_pts.insert(0, prev);
            } else {
                break;
            }
        }

        loops.push(loop_pts);
    }

    loops
}

fn draw_contours(
    chart: &mut ChartContext<
        SVGBackend,
        Cartesian2d<impl Ranged<ValueType = f64>, impl Ranged<ValueType = f64>>,
    >,
    level: f64,
    levels: &Vec<Vec<f64>>,
    edges: &Vec<Vec<(f64, f64)>>,
    grid_size: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let levels = cell_to_vertex(levels, grid_size);
    let segments = marching_squares(&levels, edges, level, grid_size);
    let loops = stitch_segments(
        &segments,
        (edges[0][0].0, edges[grid_size][grid_size].0),
        (edges[0][0].1, edges[grid_size][grid_size].1),
    );

    for loop_pts in loops {
        if loop_pts.len() < 3 {
            continue;
        }

        let mut pts: Vec<(f64, f64)> = loop_pts.iter().map(|p| (p.0, p.1)).collect();
        // pts.push(pts[0]); // close loop

        chart.draw_series(std::iter::once(PathElement::new(pts, &BLACK)))?;
    }

    Ok(())
}

fn plot_correlation(
    area: &DrawingArea<SVGBackend, Shift>,
    x_data: &[f64],
    y_data: &[f64],
    font: FontDesc<'static>,
) -> Result<(), Box<dyn std::error::Error>> {
    // Calculate Pearson correlation
    let n = x_data.len() as f64;
    let x_mean: f64 = x_data.iter().sum::<f64>() / n;
    let y_mean: f64 = y_data.iter().sum::<f64>() / n;

    let covariance: f64 = x_data
        .iter()
        .zip(y_data.iter())
        .map(|(&x, &y)| (x - x_mean) * (y - y_mean))
        .sum::<f64>()
        / n;

    let x_var: f64 = x_data.iter().map(|&x| (x - x_mean).powi(2)).sum::<f64>() / n;
    let y_var: f64 = y_data.iter().map(|&y| (y - y_mean).powi(2)).sum::<f64>() / n;

    let correlation = covariance / (x_var.sqrt() * y_var.sqrt());

    // Display correlation coefficient
    let text = format!("ρ = {:.3}", correlation);

    let mut chart = ChartBuilder::on(area)
        .margin_left(GAP)
        .margin_bottom(GAP)
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;

    let x_label_fraction = X_LABEL_HEIGHT as f64 / (area.dim_in_pixel().0 as f64);
    let y_label_fraction = Y_LABEL_WIDTH as f64 / (area.dim_in_pixel().0 as f64);
    let subplot_width = 1.0 - x_label_fraction;
    let subplot_height = 1.0 - y_label_fraction;
    let center_point = (
        0.5 * subplot_width + x_label_fraction as f64,
        0.5 * subplot_height + y_label_fraction as f64,
    );

    chart.draw_series(std::iter::once(Text::new(text, center_point, font)))?;

    // Color code by correlation strength
    let color = if correlation.abs() > 0.7 {
        RED.mix(0.3).filled()
    } else if correlation.abs() > 0.3 {
        YELLOW.mix(0.3).filled()
    } else {
        GREEN.mix(0.3).filled()
    };

    chart.draw_series(std::iter::once(Rectangle::new(
        [(0.0, 0.0), (1.0, 1.0)],
        color,
    )))?;

    Ok(())
}
