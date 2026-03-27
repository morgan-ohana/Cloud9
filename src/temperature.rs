use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A simple linear interpolator for the RELHIC T(nH) relation.
/// Stores log10(nH) and log10(T) from the digitized Benítez-Llambay+2017 Fig.5 curve.
pub struct TnRelation {
    log_nh: Vec<f64>, // log10(nH / cm^-3), sorted ascending
    log_t: Vec<f64>,  // log10(T / K)
}

impl TnRelation {
    /// Load from the CSV produced by the digitizer.
    /// Expected columns: log10_nH_cm3, log10_T_K, nH_cm3, T_K
    pub fn from_csv<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut log_nh = Vec::new();
        let mut log_t = Vec::new();

        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 {
                continue;
            } // skip header
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            let cols: Vec<&str> = line.split(',').collect();
            if cols.len() < 2 {
                return Err(format!("Line {}: expected at least 2 columns", i + 1).into());
            }
            log_nh.push(cols[0].trim().parse::<f64>()?);
            log_t.push(cols[1].trim().parse::<f64>()?);
        }

        if log_nh.is_empty() {
            return Err("CSV contained no data rows".into());
        }

        Ok(Self { log_nh, log_t })
    }

    /// Interpolate log10(T) at a given log10(nH).
    pub fn log_t_at_log_nh(&self, log_nh: f64) -> f64 {
        let n = self.log_nh.len();

        // T ~ nH^alpha, so a linear extrapolation is physically motivated.
        if log_nh <= self.log_nh[0] {
            let slope = (self.log_t[10] - self.log_t[0]) / (self.log_nh[10] - self.log_nh[0]);
            return self.log_t[0] + slope * (log_nh - self.log_nh[0]);
        }

        // Above the table: clamp — the high-density region in equilibrium so constant T expected
        if log_nh >= self.log_nh[n - 1] {
            return self.log_t[n - 1];
        }

        // Binary search for the bracketing interval
        let idx = self
            .log_nh
            .partition_point(|&x| x < log_nh)
            .saturating_sub(1);

        let x0 = self.log_nh[idx];
        let x1 = self.log_nh[idx + 1];
        let y0 = self.log_t[idx];
        let y1 = self.log_t[idx + 1];

        // Linear interpolation
        let t = (log_nh - x0) / (x1 - x0);
        y0 + t * (y1 - y0)
    }

    /// Convenience: give nH in cm^-3, get T in Kelvin.
    pub fn temperature_k(&self, nh_cm3: f64) -> f64 {
        let log_nh = nh_cm3.log10();
        let log_t = self.log_t_at_log_nh(log_nh);
        10f64.powf(log_t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plotting::plot_functions;
    use plotters::style::FontDesc;
    use plotters::style::{FontFamily, FontStyle};

    #[test]
    fn test_plot_relhic_tn_relation() -> Result<(), Box<dyn std::error::Error>> {
        let relation = TnRelation::from_csv("relhic_Tn_relation.csv")?;

        // Build points (nH, T)
        let data: Vec<(f64, f64)> = {
            let mut data = Vec::new();
            for i in 0..relation.log_nh.len() {
                data.push((
                    10_f64.powf(relation.log_nh[i]),
                    10_f64.powf(relation.log_t[i]),
                ))
            }
            data
        };

        // Build x and y range (nH & T)
        let mut x_points: Vec<f64> = Vec::new();
        let mut y_points: Vec<f64> = Vec::new();

        let (x_min, x_max): (f64, f64) = (1e-8, 1e1);

        for i in 0..1000 {
            let progress_frac = (i as f64) / 999.0;
            x_points.push((x_min.ln() + progress_frac * (x_max.ln() - x_min.ln())).exp());
            y_points.push(relation.temperature_k(*x_points.last().unwrap()));
        }

        let font = FontDesc::new(FontFamily::SansSerif, 16.0, FontStyle::Normal);

        plot_functions(
            &x_points,
            &vec![y_points],
            "relhic_tn_relation.png",
            "RELHIC T–nH Relation (Benítez-Llambay+2017 Fig.5)",
            "log₁₀(nH / cm⁻³)",
            "log₁₀(T / K)",
            vec![Some("T(nH) composite curve".to_string())],
            font,
            vec![true],  // dashed
            Some(&data), // no scatter data
            None,        // no y error bars
        )?;

        Ok(())
    }
}
