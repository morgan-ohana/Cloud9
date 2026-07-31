use plotters::coord::ranged1d::ValueFormatter;
use plotters::prelude::*;

pub fn linear_ticks(min: f64, max: f64, n: usize) -> Vec<f64> {
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

pub fn log_ticks(min: f64, max: f64, level: i32) -> Vec<f64> {
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

pub fn draw_ticks_top_and_right<
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
