use plotters::prelude::*;
use plotters::coord::Shift;

pub fn plot_function(
    x_points: &Vec<f64>,
    y_points: &Vec<f64>,
    filename: &str,
    title: &str,
    xlabel: &str,
    ylabel: &str,
    legend: Option<String>,
    data: Option<&Vec<(f64, f64)>>
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
        },
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
        Some(data_points) => 0.9*data_points[0].0..1.1*data_points.last().unwrap().0,
        None => x_points[0]..x_points[x_points.len() - 1]
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
        chart.draw_series(LineSeries::new(plot_profile, &BLUE))?.label(legend_text);
        // Configure and draw legend
        chart.configure_series_labels()
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw()?;
    } else {
        chart.draw_series(LineSeries::new(plot_profile, &BLUE))?;
    }

    if let Some(data_points) = data {
        chart.draw_series(data_points.iter().map(|point| {
            Circle::new(*point, 5, &RED)
        })).unwrap();
    }

    root.present()?;
    //println!("Plot saved as {}", filename);
    Ok(())
}

const LABEL_WIDTH: u32 = 30;
pub fn create_corner_plot(
    chain: &[[f64; 4]],
    param_names: &[&str; 4],
    output_path: &str,
    burn_in_fraction: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    // Remove burn-in
    let burn_in = (chain.len() as f64 * burn_in_fraction) as usize;
    let chain = &chain[burn_in..];
    
    // Extract parameter columns
    let mut params: Vec<Vec<f64>> = vec![Vec::new(); 4];
    for point in chain {
        for i in 0..4 {
            params[i].push(point[i]);
        }
    }
    
    // Create plot area
    let root = BitMapBackend::new(output_path, (1600, 1600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Split into 4x4 subplots
    let sub_areas = root
        .margin(5, 50, 50, 5)
        .split_evenly((4, 4));
    
    // Plot each cell
    for row in 0..4 {
        for col in 0..4 {
            let idx = row * 4 + col;
            let drawing_area = &sub_areas[idx];
            
            if row == col {
                // Diagonal: Histogram
                plot_histogram(drawing_area, &params[row], param_names[row], (row, col))?;
            } else if row > col {
                // Lower triangle: 2D scatter/density
                plot_2d_scatter(drawing_area, &params[col], &params[row], param_names, (row, col))?;
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
    (row, col): (usize, usize)
) -> Result<(), Box<dyn std::error::Error>> {
    // Calculate bins
    let min = data.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let max = data.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    let n_bins = 50;
    let bin_width = (max - min) / n_bins as f64;
    let mut bins = vec![0; n_bins];
    
    for &value in data {
        let bin_idx = ((value - min) / bin_width).floor() as usize;
        let bin_idx = bin_idx.min(n_bins - 1);
        bins[bin_idx] += 1;
    }
    
    let max_count = *bins.iter().max().unwrap() as f64;
    
    // Create chart for histogram
    let mut chart_builder = ChartBuilder::on(area);
    
    if col == 0 {
        chart_builder.y_label_area_size(LABEL_WIDTH);
    } else {
        chart_builder.margin_left(LABEL_WIDTH);
    }

    if row == 3 {
        chart_builder.x_label_area_size(LABEL_WIDTH);
    } else {
        chart_builder.margin_bottom(LABEL_WIDTH);
    }

    let mut chart = chart_builder
        //.caption(param_name, ("sans-serif", 15).into_font())
        .build_cartesian_2d(min..max, 0.0..max_count * 1.1)?;
    
    chart.configure_mesh()
        .x_desc(param_name) // X-axis label
        .y_desc("Counts") // Y-axis label
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
        }).draw()?;
    
    // Plot histogram bars
    for i in 0..n_bins {
        let x_start = min + i as f64 * bin_width;
        let x_end = x_start + bin_width;
        let count = bins[i] as f64;
        
        chart.draw_series(std::iter::once(Rectangle::new(
            [(x_start, 0.0), (x_end, count)],
            BLUE.mix(0.5).filled(),
        )))?;
    }
    
    // Add KDE curve
    plot_kde(area, data, min, max)?;
    
    Ok(())
}

fn plot_2d_scatter(
    area: &DrawingArea<BitMapBackend, Shift>,
    x_data: &[f64],
    y_data: &[f64],
    param_names: &[&str; 4],
    (row, col): (usize, usize)
) -> Result<(), Box<dyn std::error::Error>> {
    // Thin the data if too many points
    let thin_factor = (x_data.len() / 5000).max(1);
    let thinned_x: Vec<f64> = x_data.iter().step_by(thin_factor).copied().collect();
    let thinned_y: Vec<f64> = y_data.iter().step_by(thin_factor).copied().collect();
    
    let x_min = thinned_x.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let x_max = thinned_x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let y_min = thinned_y.iter().fold(f64::INFINITY, |a, &b| a.min(b));
    let y_max = thinned_y.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    let mut chart_builder = ChartBuilder::on(area);
    
    if col == 0 {
        chart_builder.y_label_area_size(LABEL_WIDTH);
    } else {
        chart_builder.margin_left(LABEL_WIDTH);
    }

    if row == 3 {
        chart_builder.x_label_area_size(LABEL_WIDTH);
    } else {
        chart_builder.margin_bottom(LABEL_WIDTH);
    }

    let mut chart = chart_builder.build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    
    let mut mesh = chart.configure_mesh();
    
    mesh.x_desc(param_names[col]) // X-axis label
        .y_desc(param_names[row]) // Y-axis label
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
        });
    
    mesh.draw()?;
    
    // Create scatter plot with transparency
    chart.draw_series(
        thinned_x.iter().zip(thinned_y.iter()).enumerate().map(|(idx,(&x, &y))| {
            let t = idx as f64 / thinned_x.len() as f64;
            
            // Color interpolation: Blue (early) -> Purple (middle) -> Red (late)
            let color = if t < 0.5 {
                // Blue to Purple
                let u = t * 2.0; // 0 to 1
                RGBColor(
                    (255.0 * u) as u8,      // R increases
                    0,                      // G stays 0
                    (255.0 * (1.0 - u)) as u8, // B decreases
                )
            } else {
                // Purple to Red
                let u = (t - 0.5) * 2.0; // 0 to 1
                RGBColor(
                    255,                    // R stays max
                    (128.0 * u) as u8,     // G increases slightly
                    (255.0 * (1.0 - u)) as u8, // B decreases
                )
            };

            Circle::new((x, y), 1, color.mix(0.1).filled())
        })
    )?;
    
    // Add contour lines for density
    plot_2d_contours(area, x_data, y_data, x_min, x_max, y_min, y_max)?;
    
    Ok(())
}

fn plot_kde(
    area: &DrawingArea<BitMapBackend, Shift>,
    data: &[f64],
    min: f64,
    max: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    // Simple KDE using Gaussian kernel
    let bandwidth = (max - min) / 50.0;
    let n_points = 200;
    
    let mut kde_points = Vec::new();
    for i in 0..n_points {
        let x = min + (max - min) * i as f64 / n_points as f64;
        
        let mut density = 0.0;
        for &point in data {
            let diff = (x - point) / bandwidth;
            density += (-0.5 * diff * diff).exp();
        }
        
        density /= (data.len() as f64 * bandwidth * (2.0 * std::f64::consts::PI).sqrt());
        kde_points.push((x, density));
    }
    
    // Scale to match histogram
    let max_density = kde_points.iter().map(|&(_, d)| d).fold(f64::NEG_INFINITY, f64::max);
    let scale_factor = {
        let mut hist_area = area.clone();
        let hist_chart = ChartBuilder::on(&mut hist_area)
            .build_cartesian_2d(min..max, 0.0..1.0)?;
        let y_scale = 1.0 / max_density;
        y_scale
    };
    
    let mut chart = ChartBuilder::on(area)
        .margin_left(LABEL_WIDTH).margin_bottom(LABEL_WIDTH)
        .build_cartesian_2d(min..max, 0.0..max_density * scale_factor)?;
    
    chart.draw_series(LineSeries::new(
        kde_points.iter().map(|&(x, d)| (x, d * scale_factor)),
        &RED,
    ))?;
    
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
) -> Result<(), Box<dyn std::error::Error>> {
    // Simple 2D density estimation
    let grid_size = 50;
    let x_bin_width = (x_max - x_min) / grid_size as f64;
    let y_bin_width = (y_max - y_min) / grid_size as f64;
    
    let mut density = vec![vec![0.0; grid_size]; grid_size];
    
    for (&x, &y) in x_data.iter().zip(y_data.iter()) {
        let x_bin = ((x - x_min) / x_bin_width).floor() as usize;
        let y_bin = ((y - y_min) / y_bin_width).floor() as usize;
        
        if x_bin < grid_size && y_bin < grid_size {
            density[x_bin][y_bin] += 1.0;
        }
    }
    
    // Normalize
    let max_density = density.iter()
        .flatten()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    
    if max_density > 0.0 {
        for row in &mut density {
            for val in row {
                *val /= max_density;
            }
        }
        
        // Draw contour lines at 10%, 30%, 50%, 70%, 90%
        let contours = [0.1, 0.3, 0.5, 0.7, 0.9];
        
        for &level in &contours {
            let mut points = Vec::new();
            
            for i in 0..grid_size - 1 {
                for j in 0..grid_size - 1 {
                    let corners = [
                        density[i][j],
                        density[i+1][j],
                        density[i][j+1],
                        density[i+1][j+1],
                    ];
                    
                    // Simple contour detection (Marching Squares)
                    if corners.iter().any(|&d| d >= level) && 
                       corners.iter().any(|&d| d < level) {
                        let x = x_min + (i as f64 + 0.5) * x_bin_width;
                        let y = y_min + (j as f64 + 0.5) * y_bin_width;
                        points.push((x, y));
                    }
                }
            }
            
            if !points.is_empty() {
                let mut chart = ChartBuilder::on(area)
                    .margin_left(LABEL_WIDTH).margin_bottom(LABEL_WIDTH)
                    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
                
                chart.draw_series(
                    AreaSeries::new(points, 0.0, BLUE.mix(level as f64 * 0.3))
                )?;
            }
        }
    }
    
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
    
    let covariance: f64 = x_data.iter().zip(y_data.iter())
        .map(|(&x, &y)| (x - x_mean) * (y - y_mean))
        .sum::<f64>() / n;
    
    let x_var: f64 = x_data.iter().map(|&x| (x - x_mean).powi(2)).sum::<f64>() / n;
    let y_var: f64 = y_data.iter().map(|&y| (y - y_mean).powi(2)).sum::<f64>() / n;
    
    let correlation = covariance / (x_var.sqrt() * y_var.sqrt());
    
    // Display correlation coefficient
    let text = format!("ρ = {:.3}", correlation);
    
    let mut chart = ChartBuilder::on(area)
        .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;
    
    chart.draw_series(std::iter::once(
        Text::new(
            text,
            (0.5, 0.5),
            ("sans-serif", 20).into_font()
        )
    ))?;
    
    // Color code by correlation strength
    let color = if correlation.abs() > 0.7 {
        RED.mix(0.3).filled()
    } else if correlation.abs() > 0.3 {
        YELLOW.mix(0.3).filled()
    } else {
        GREEN.mix(0.3).filled()
    };
    
    chart.draw_series(std::iter::once(Rectangle::new([(0.0, 0.0), (1.0, 1.0)], color)))?;
    
    Ok(())
}