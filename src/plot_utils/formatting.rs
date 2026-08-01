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
