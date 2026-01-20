use plotters::prelude::*;

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
