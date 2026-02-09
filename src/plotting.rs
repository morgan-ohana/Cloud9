use plotters::coord::Shift;
use plotters::coord::ranged1d::ValueFormatter;
use plotters::prelude::*;

pub fn plot_function(
    x_points: &Vec<f64>,
    y_points: &Vec<f64>,
    filename: &str,
    title: &str,
    xlabel: &str,
    ylabel: &str,
    legend: Option<String>,
    data: Option<&Vec<(f64, f64)>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;

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
        }
        None => {
            for i in 0..y_points.len() {
                if y_points[i] > y_max {
                    y_max = y_points[i]
                }
                if y_points[i] < y_min {
                    y_min = y_points[i]
                }
            }
        }
    }

    //println!("y_min = {:.3}, y_max={:.3}", y_min, y_max);

    let x_range = match data {
        Some(data_points) => 0.9 * data_points[0].0..1.1 * data_points.last().unwrap().0,
        None => x_points[0]..x_points[x_points.len() - 1],
    };

    let x_range = x_range.log_scale();

    let y_range = (y_min + 1e-4)
        * match y_min.signum() {
            1.0 => 0.9,
            -1.0 => 1.1,
            _ => panic!("number has no sign, is probably NaN"),
        }
        ..y_max
            * match y_max.signum() {
                1.0 => 1.1,
                -1.0 => 0.9,
                _ => panic!("number has no sign, is probably NaN"),
            };

    let y_range = y_range.log_scale();

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(x_range, y_range)?;

    chart
        .configure_mesh()
        .x_desc(xlabel) // X-axis label
        .y_desc(ylabel) // Y-axis label
        .x_label_formatter(&|x| {
            if x.abs() >= 1000.0 || x.abs() <= 0.1 {
                format!("{:.1e}", x)
            } else {
                format!("{:.1}", x)
            }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    let plot_profile: Vec<(f64, f64)> = (0..x_points.len())
        .map(|i| (x_points[i], y_points[i]))
        .collect();

    if let Some(legend_text) = legend {
        chart
            .draw_series(LineSeries::new(plot_profile, &BLUE))?
            .label(legend_text);
        // Configure and draw legend
        chart
            .configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;
    } else {
        chart.draw_series(LineSeries::new(plot_profile, &BLUE))?;
    }

    if let Some(data_points) = data {
        chart
            .draw_series(data_points.iter().map(|point| Circle::new(*point, 5, &RED)))
            .unwrap();
    }

    root.present()?;
    //println!("Plot saved as {}", filename);
    Ok(())
}

pub fn create_chain_trace_plots(
    chains: &[([f64; 4], Vec<[f64; 4]>, Vec<f64>)],
    burn_in: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    for (chain_id, (_, chain, _)) in chains.iter().enumerate() {
        let filename = format!("trace_plots/chain_{}_trace.png", chain_id);
        let root = BitMapBackend::new(&filename, (1200, 800)).into_drawing_area();
        root.fill(&WHITE)?;

        // Create subplots for each parameter
        let sub_areas = root.split_evenly((2, 2));
        let param_names = ["m200_0", "c200_0", "tau", "rho_c"];

        for (param_idx, area) in sub_areas.into_iter().enumerate() {
            if param_idx >= 4 {
                break;
            }

            let data: Vec<f64> = chain.iter().map(|p| p[param_idx]).collect();

            // Find reasonable y-range
            let mut sorted: Vec<f64> = data.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let p05 = sorted[(sorted.len() as f64 * 0.05) as usize];
            let p95 = sorted[(sorted.len() as f64 * 0.95) as usize];
            let y_range = p05..p95;

            let mut chart = ChartBuilder::on(&area)
                .caption(
                    format!("Chain {}: {}", chain_id, param_names[param_idx]),
                    ("sans-serif", 16),
                )
                .margin(10)
                .build_cartesian_2d(0..chain.len(), y_range)?;

            chart
                .configure_mesh()
                .x_desc("Step")
                .y_desc("Value")
                .label_style(("sans-serif", 10))
                .draw()?;

            // Plot trace
            chart.draw_series(LineSeries::new(
                data.iter().enumerate().map(|(i, &v)| (i, v)),
                &BLUE,
            ))?;

            // Mark burn-in cutoff
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(burn_in, p05), (burn_in, p95)],
                &RED,
            )))?;
        }

        root.present()?;
        println!("Saved trace plot: {}", filename);
    }

    Ok(())
}

const LABEL_WIDTH: u32 = 30;
const Y_LABEL_PAD: u32 = 15;
const LOG_SCALE: [bool; 4] = [true, false, false, true];
pub fn create_corner_plot(
    chain: &[[f64; 4]],
    marked_points: &[&[f64; 4]],
    param_names: &[&str; 4],
    output_path: &str,
    burn_in_fraction: f64,
    inwards: bool,
    bounds: &[[f64; 2]; 4],
) -> Result<(), Box<dyn std::error::Error>> {
    // Remove burn-in
    let burn_in = (chain.len() as f64 * burn_in_fraction) as usize;
    let chain = &chain[burn_in..];

    let num_params = match inwards {
        true => 3,
        false => 4,
    };

    // Extract parameter columns
    let mut params: Vec<Vec<f64>> = vec![Vec::new(); num_params];
    for point in chain {
        for i in 0..num_params {
            params[i].push(point[i]);
        }
    }

    // Create plot area
    let root = BitMapBackend::new(output_path, (1600, 1600)).into_drawing_area();
    root.fill(&WHITE)?;

    let plot_width = 1600 - LABEL_WIDTH - 5;

    let x_break_points = [plot_width/4 + Y_LABEL_PAD, plot_width/2 + Y_LABEL_PAD, 3*plot_width/4 + Y_LABEL_PAD];
    let y_break_points = [plot_width/4, plot_width/2, 3*plot_width/4];

    // Split into 4x4 subplots
    let sub_areas = root
        .margin(5, 5 + LABEL_WIDTH, 5 + LABEL_WIDTH - Y_LABEL_PAD, 5)
        .split_by_breakpoints(x_break_points, y_break_points);
        //.split_evenly((num_params, num_params));

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
                    (row, col),
                    &bounds[row],
                    &marked_values,
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
                    (row, col),
                    &bounds,
                    &marked_2d,
                )?;
            } else {
                // Upper triangle: Correlation/contour or leave empty
                plot_correlation(drawing_area, &params[col], &params[row])?;
            }
        }
    }

    root.present()?;
    println!("Corner plot saved to: {}", output_path);
    Ok(())
}

fn plot_histogram(
    area: &DrawingArea<BitMapBackend, Shift>,
    data: &[f64],
    param_name: &str,
    (row, col): (usize, usize),
    bounds: &[f64; 2],
    marked_values: &[f64],
) -> Result<(), Box<dyn std::error::Error>> {
    // Calculate bins
    let min = bounds[0];
    let max = bounds[1];
    // let min = data.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    // let max = data.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    let n_bins = 50;
    let edges = {
        match LOG_SCALE[row] {
            true => {
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

    let spacing = match LOG_SCALE[row] {
        true => (max.ln() - min.ln()) / (n_bins as f64),
        false => (max - min) / (n_bins as f64),
    };
    let mut bins: Vec<u32> = vec![0; n_bins];

    for &value in data {
        if value >= max || value <= min {
            continue;
        }
        let bin_idx = match LOG_SCALE[row] {
            true => ((value.ln() - min.ln()) / spacing).floor() as usize,
            false => ((value - min) / spacing).floor() as usize,
        };
        bins[bin_idx] += 1;
    }

    let max_count = *bins.iter().max().unwrap() as f64;

    // Create chart for histogram
    let mut chart_builder = ChartBuilder::on(area);

    if col == 0 {
        chart_builder.y_label_area_size(LABEL_WIDTH + Y_LABEL_PAD);
    } else {
        chart_builder.margin_left(LABEL_WIDTH);
    }

    if row == 3 {
        chart_builder.x_label_area_size(LABEL_WIDTH);
    } else {
        chart_builder.margin_bottom(LABEL_WIDTH);
    }

    match LOG_SCALE[row] {
        true => {
            let mut chart = chart_builder
                //.caption(param_name, ("sans-serif", 15).into_font())
                .build_cartesian_2d((min..max).log_scale(), 0.0..max_count * 1.1)?;

            draw_hist_content(&mut chart, param_name, &bins, &edges)?;

            // Draw marked values with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &marked_value) in marked_values.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(PathElement::new(
                    vec![(marked_value, 0.0), (marked_value, max_count * 1.1)],
                    color,
                )))?;
            }
        }
        false => {
            let mut chart = chart_builder
                //.caption(param_name, ("sans-serif", 15).into_font())
                .build_cartesian_2d(min..max, 0.0..max_count * 1.1)?;

            draw_hist_content(&mut chart, param_name, &bins, &edges)?;

            // Draw marked values with different colors
            let colors = [&GREEN, &RED, &BLUE, &MAGENTA];
            for (i, &marked_value) in marked_values.iter().enumerate() {
                let color = colors[i % colors.len()];
                chart.draw_series(std::iter::once(PathElement::new(
                    vec![(marked_value, 0.0), (marked_value, max_count * 1.1)],
                    color,
                )))?;
            }
        }
    };

    // Add KDE curve
    plot_kde(area, data, min, max, (row, col), LOG_SCALE[row])?;

    Ok(())
}

fn draw_hist_content<
    X: Ranged<ValueType = f64> + ValueFormatter<f64>,
    Y: Ranged<ValueType = f64> + ValueFormatter<f64>,
>(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<X, Y>>,
    param_name: &str,
    bins: &Vec<u32>,
    edges: &Vec<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    chart
        .configure_mesh()
        .x_desc(param_name) // X-axis label
        .y_desc("Counts") // Y-axis label
        .x_label_formatter(&|x| {
            if x.abs() >= 1000.0 || x.abs() <= 0.01 {
                format!("{:.1e}", x)
            } else {
                format!("{:.1}", x)
            }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.01 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    // Plot histogram bars
    for i in 0..bins.len() {
        let count = bins[i] as f64;

        chart.draw_series(std::iter::once(Rectangle::new(
            [(edges[i], 0.0), (edges[i + 1], count)],
            BLUE.mix(0.5).filled(),
        )))?;
    }

    Ok(())
}

fn plot_2d_scatter(
    area: &DrawingArea<BitMapBackend, Shift>,
    x_data: &[f64],
    y_data: &[f64],
    param_names: &[&str; 4],
    (row, col): (usize, usize),
    bounds: &[[f64; 2]; 4],
    marked_points_2d: &[(f64, f64)],
) -> Result<(), Box<dyn std::error::Error>> {
    // Thin the data if too many points
    let thin_factor = (x_data.len() / 5000).max(1);
    let thinned_x: Vec<f64> = x_data.iter().step_by(thin_factor).copied().collect();
    let thinned_y: Vec<f64> = y_data.iter().step_by(thin_factor).copied().collect();

    let [x_min, x_max] = bounds[col];
    let [y_min, y_max] = bounds[row];

    let mut chart_builder = ChartBuilder::on(area);

    if col == 0 {
        chart_builder.y_label_area_size(LABEL_WIDTH + Y_LABEL_PAD);
    } else {
        chart_builder.margin_left(LABEL_WIDTH);
    }

    if row == 3 {
        chart_builder.x_label_area_size(LABEL_WIDTH);
    } else {
        chart_builder.margin_bottom(LABEL_WIDTH);
    }

    match (LOG_SCALE[col], LOG_SCALE[row]) {
        (true, true) => {
            let mut chart = chart_builder
                .build_cartesian_2d((x_min..x_max).log_scale(), (y_min..y_max).log_scale())?;

            draw_scatter_content(&mut chart, &thinned_x, &thinned_y, param_names, (row, col))?;

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

            draw_scatter_content(&mut chart, &thinned_x, &thinned_y, param_names, (row, col))?;

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

            draw_scatter_content(&mut chart, &thinned_x, &thinned_y, param_names, (row, col))?;

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

            draw_scatter_content(&mut chart, &thinned_x, &thinned_y, param_names, (row, col))?;

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
    plot_2d_contours(area, x_data, y_data, x_min, x_max, y_min, y_max, (row, col))?;

    Ok(())
}

fn draw_scatter_content<
    X: Ranged<ValueType = f64> + ValueFormatter<f64>,
    Y: Ranged<ValueType = f64> + ValueFormatter<f64>,
>(
    chart: &mut ChartContext<BitMapBackend, Cartesian2d<X, Y>>,
    thinned_x: &Vec<f64>,
    thinned_y: &Vec<f64>,
    param_names: &[&str; 4],
    (row, col): (usize, usize),
) -> Result<(), Box<dyn std::error::Error>> {
    let mut mesh = chart.configure_mesh();

    mesh.x_desc(param_names[col]) // X-axis label
        .y_desc(param_names[row]) // Y-axis label
        .x_label_formatter(&|x| {
            if x.abs() >= 1000.0 || x.abs() <= 0.01 {
                format!("{:.1e}", x)
            } else {
                format!("{:.1}", x)
            }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.01 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        });

    mesh.draw()?;

    // Create scatter plot with transparency
    chart.draw_series(thinned_x.iter().zip(thinned_y.iter()).enumerate().map(
        |(idx, (&x, &y))| {
            let t = idx as f64 / thinned_x.len() as f64;

            // Color interpolation: Blue (early) -> Purple (middle) -> Red (late)
            let color = if t < 0.5 {
                // Blue to Purple
                let u = t * 2.0; // 0 to 1
                RGBColor(
                    (255.0 * u) as u8,         // R increases
                    0,                         // G stays 0
                    (255.0 * (1.0 - u)) as u8, // B decreases
                )
            } else {
                // Purple to Red
                let u = (t - 0.5) * 2.0; // 0 to 1
                RGBColor(
                    255,                       // R stays max
                    (128.0 * u) as u8,         // G increases slightly
                    (255.0 * (1.0 - u)) as u8, // B decreases
                )
            };

            Circle::new((x, y), 1, color.mix(0.1).filled())
        },
    ))?;

    Ok(())
}

fn plot_kde(
    area: &DrawingArea<BitMapBackend, Shift>,
    data: &[f64],
    min: f64,
    max: f64,
    (row, col): (usize, usize),
    log_scale: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Simple KDE using Gaussian kernel
    let bandwidth = match log_scale {
        true => (max.ln() - min.ln()) / 50.0,
        false => (max - min) / 50.0,
    };
    let n_points = 200;

    let mut kde_points = Vec::new();
    for i in 0..n_points {
        let x = match log_scale {
            true => (min.ln() + (max.ln() - min.ln()) * i as f64 / n_points as f64).exp(),
            false => min + (max - min) * i as f64 / n_points as f64,
        };

        let mut density = 0.0;
        for &point in data {
            let diff = match log_scale {
                true => (x.ln() - point.ln()) / bandwidth,
                false => (x - point) / bandwidth,
            };
            density += (-0.5 * diff * diff).exp();
        }

        density /= data.len() as f64 * bandwidth * (2.0 * std::f64::consts::PI).sqrt();
        kde_points.push((x, density));
    }

    // Scale to match histogram
    let max_density = kde_points
        .iter()
        .map(|&(_, d)| d)
        .fold(f64::NEG_INFINITY, f64::max);
    let scale_factor = 1.0 / (1.1 * max_density);

    let y_margin_width = match col {
        0 => LABEL_WIDTH + Y_LABEL_PAD,
        _ => LABEL_WIDTH
    };

    match log_scale {
        true => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(LABEL_WIDTH)
                .build_cartesian_2d((min..max).log_scale(), 0.0..1.0)?;

            chart.draw_series(LineSeries::new(
                kde_points.iter().map(|&(x, d)| (x, d * scale_factor)),
                &RED,
            ))?;
        }
        false => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(LABEL_WIDTH)
                .build_cartesian_2d(min..max, 0.0..1.0)?;

            chart.draw_series(LineSeries::new(
                kde_points.iter().map(|&(x, d)| (x, d * scale_factor)),
                &RED,
            ))?;
        }
    }

    Ok(())
}

fn plot_2d_contours(
    area: &DrawingArea<BitMapBackend, Shift>,
    x_data: &[f64],
    y_data: &[f64],
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    (row, col): (usize, usize),
) -> Result<(), Box<dyn std::error::Error>> {
    // Simple 2D density estimation
    let grid_size = 50;

    let mut density = vec![vec![0.0; grid_size]; grid_size];
    let mut edges = vec![vec![(0.0, 0.0); grid_size + 1]; grid_size + 1];
    for i in 0..=grid_size {
        for j in 0..=grid_size {
            let x_edge = match LOG_SCALE[col] {
                true => {
                    (x_min.ln() + (x_max.ln() - x_min.ln()) * i as f64 / grid_size as f64).exp()
                }
                false => x_min + (x_max - x_min) * i as f64 / grid_size as f64,
            };
            let y_edge = match LOG_SCALE[row] {
                true => {
                    (y_min.ln() + (y_max.ln() - y_min.ln()) * j as f64 / grid_size as f64).exp()
                }
                false => y_min + (y_max - y_min) * j as f64 / grid_size as f64,
            };
            edges[i][j] = (x_edge, y_edge);
        }
    }

    let mut count = 0;
    for (&x, &y) in x_data.iter().zip(y_data.iter()) {
        if x > x_max || x < x_min || y > y_max || y < y_min {
            continue;
        }

        let x_bin = match LOG_SCALE[col] {
            true => ((x.ln() - x_min.ln()) / (x_max.ln() - x_min.ln()) * grid_size as f64).floor()
                as usize,
            false => ((x - x_min) / (x_max - x_min) * grid_size as f64).floor() as usize,
        };
        let y_bin = match LOG_SCALE[row] {
            true => ((y.ln() - y_min.ln()) / (y_max.ln() - y_min.ln()) * grid_size as f64).floor()
                as usize,
            false => ((y - y_min) / (y_max - y_min) * grid_size as f64).floor() as usize,
        };

        count += 1;
        density[x_bin][y_bin] += 1.0;
    }

    //Normalize
    for i in 0..grid_size {
        for j in 0..grid_size {
            density[i][j] /= count as f64
        }
    }

    let y_margin_width = match col {
        0 => LABEL_WIDTH + Y_LABEL_PAD,
        _ => LABEL_WIDTH
    };

    match (LOG_SCALE[col], LOG_SCALE[row]) {
        (true, true) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(LABEL_WIDTH)
                .build_cartesian_2d((x_min..x_max).log_scale(), (y_min..y_max).log_scale())?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
        (true, false) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(LABEL_WIDTH)
                .build_cartesian_2d((x_min..x_max).log_scale(), y_min..y_max)?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
        (false, true) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(LABEL_WIDTH)
                .build_cartesian_2d(x_min..x_max, (y_min..y_max).log_scale())?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
        (false, false) => {
            let mut chart = ChartBuilder::on(area)
                .margin_left(y_margin_width)
                .margin_bottom(LABEL_WIDTH)
                .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

            draw_contour_content(&mut chart, &density, &edges, grid_size)?;
        }
    }

    Ok(())
}

fn draw_contour_content(
    chart: &mut ChartContext<
        BitMapBackend,
        Cartesian2d<impl Ranged<ValueType = f64>, impl Ranged<ValueType = f64>>,
    >,
    density: &Vec<Vec<f64>>,
    edges: &Vec<Vec<(f64, f64)>>,
    grid_size: usize,
) -> Result<(), Box<dyn std::error::Error>> {

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
    let mut cumulative_probability = 0.0;
    let mut levels = vec![vec![0; grid_size]; grid_size];
    
    for mut cell in cells {
        cumulative_probability += cell.0;
        
        let color = if cumulative_probability <= CONTOUR_LEVELS[0] {
            levels[cell.3.0][cell.3.1] = 0;
            RGBColor(255, 100, 100) // Red
        } else if cumulative_probability <= CONTOUR_LEVELS[1] {
            levels[cell.3.0][cell.3.1] = 1;
            RGBColor(255, 200, 100) // Orange
        } else if cumulative_probability <= CONTOUR_LEVELS[2] {
            levels[cell.3.0][cell.3.1] = 2;
            RGBColor(100, 255, 200) // Cyan
        } else if cumulative_probability <= CONTOUR_LEVELS[3] {
            levels[cell.3.0][cell.3.1] = 3;
            RGBColor(100, 200, 255) // Light Blue
        } else if cumulative_probability <= CONTOUR_LEVELS[4] {
            levels[cell.3.0][cell.3.1] = 4;
            RGBColor(100, 100, 255) // Blue
        } else {
            levels[cell.3.0][cell.3.1] = 5;
            continue;
        };
    
        chart.draw_series(std::iter::once(Rectangle::new(
            [cell.1, cell.2],
            color.mix(0.3).filled(),
        )))?;
    }

    // Draw contours
    draw_contour_around_level(chart, 1, &levels, edges, grid_size)?;
    draw_contour_around_level(chart, 3, &levels, edges, grid_size)?;

    Ok(())
}

fn draw_contour_around_level(
    chart: &mut ChartContext<
        BitMapBackend,
        Cartesian2d<impl Ranged<ValueType = f64>, impl Ranged<ValueType = f64>>,
    >,
    level: u32,
    levels: &Vec<Vec<u32>>,
    edges: &Vec<Vec<(f64, f64)>>,
    grid_size: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut points = Vec::new();

    for i in 0..grid_size {
        for j in 0..grid_size {
            for a in i.saturating_sub(1)..=(i+1).min(grid_size - 1) {
                for b in j.saturating_sub(1)..=(j+1).min(grid_size - 1) {
                    if levels[i][j] != level {
                        continue;
                    }
                    
                    // contour edge in plot
                    if levels[a][b] > level {
                        // edge i lies between i and i-1 so for two neighbors the greater is the idx of the edge between them
                        points.push(edges[i.max(a)][j.max(b)]);
                    }

                    // plot edge
                    if i == 0 {
                        points.push(edges[0][j]);
                        points.push(edges[0][j+1]);
                    } else if i == grid_size - 1 {
                        points.push(edges[grid_size][j]);
                        points.push(edges[grid_size][j+1]);
                    }

                    if j == 0 {
                        points.push(edges[i][0]);
                        points.push(edges[i+1][0]);
                    } else if j == grid_size - 1 {
                        points.push(edges[i][grid_size]);
                        points.push(edges[i+1][grid_size]);
                    }
                }
            }
        }
    }

    // Calculate centroid
    let (sum_x, sum_y) = points.iter()
        .fold((0.0, 0.0), |(sx, sy), &(x, y)| (sx + x, sy + y));
    let centroid = (sum_x / points.len() as f64, sum_y / points.len() as f64);
    
    // Sort by angle relative to centroid
    points.sort_by(|&(x1, y1), &(x2, y2)| {
        let angle1 = (y1 - centroid.1).atan2(x1 - centroid.0);
        let angle2 = (y2 - centroid.1).atan2(x2 - centroid.0);
        angle1.partial_cmp(&angle2).unwrap_or(std::cmp::Ordering::Equal)
    });

    chart.draw_series(std::iter::once(PathElement::new(
        points,
        &BLACK,
    )))?;

    Ok(())
}

fn plot_correlation(
    area: &DrawingArea<BitMapBackend, Shift>,
    x_data: &[f64],
    y_data: &[f64],
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
        .margin_left(LABEL_WIDTH)
        .margin_bottom(LABEL_WIDTH)
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;

    let label_fraction = LABEL_WIDTH as f64 / (area.dim_in_pixel().0 as f64);
    let subplot_width = 1.0 - label_fraction;
    let center_point = (
        0.5 * subplot_width + label_fraction as f64,
        0.5 * subplot_width + label_fraction as f64,
    );

    chart.draw_series(std::iter::once(Text::new(
        text,
        center_point,
        ("sans-serif", 20).into_font(),
    )))?;

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
