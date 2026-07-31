pub fn window_series(
    pts: &[(f64, f64)],
    x_range: (f64, f64),
    y_range: (f64, f64),
) -> Vec<(f64, f64)> {
    let (x_min, x_max) = x_range;
    let (y_min, y_max) = y_range;
    pts.iter()
        .copied()
        .filter(|&(x, y)| {
            x.is_finite()
                && y.is_finite()
                && y > 0.0 // ln(0) = -inf on a log axis -> garbage pixel coords
                && x >= x_min && x <= x_max
                && y >= y_min && y <= y_max
        })
        .collect()
}
