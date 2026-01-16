use plotters::prelude::*;

pub fn plot_function(
    x_points: &Vec<f64>,
    y_points: &Vec<f64>,
    filename: &str,
    title: &str,
    xlabel: &str,
    ylabel: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;

    for i in 0..y_points.len() {
        if y_points[i] > y_max {
            y_max = y_points[i]
        }
        if y_points[i] < y_min {
            y_min = y_points[i]
        }
    }

    println!("y_min = {:.3}, y_max={:.3}", y_min, y_max);

    let x_range = (x_points[0]..x_points[x_points.len() - 1]).log_scale();

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

    chart.draw_series(LineSeries::new(plot_profile, &BLUE))?;

    root.present()?;
    println!("Plot saved as {}", filename);
    Ok(())
}
