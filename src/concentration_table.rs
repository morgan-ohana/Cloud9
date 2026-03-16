/// Mass-concentration lookup table generated from Diemer (2019) via colossus.
///
/// The binary file layout (little-endian) is:
///   [0..4)         u32     N  (number of sample points)
///   [4..4+8N)      f64[N]  log10(M / [Msun/h]) – uniformly spaced, ascending
///   [4+8N..4+16N)  f64[N]  concentration c
///
/// Interpolation is performed in log10(M) vs log10(c) space, which is
/// appropriate because the relation is close to a power law.
use std::fs;
use std::path::Path;

/// A loaded and ready-to-query mass–concentration table.
pub struct MassConcentrationTable {
    /// log10(M / [Msun/h]) sample points (length N, ascending)
    log_mass: Vec<f64>,
    /// log10(concentration) at each sample point (length N)
    log_conc: Vec<f64>,
}

impl MassConcentrationTable {
    /// Load the lookup table produced by `gen_concentration_table.py`.
    ///
    /// # Errors
    /// Returns an error string if the file cannot be read or is malformed.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let bytes = fs::read(path).map_err(|e| format!("IO error: {e}"))?;

        // ── Parse header ────────────────────────────────────────────────────
        if bytes.len() < 4 {
            return Err("File too short to contain header".into());
        }
        let n = u32::from_le_bytes(bytes[0..4].try_into().unwrap()) as usize;

        let expected = 4 + 16 * n; // 4-byte header + two f64 arrays of length n
        if bytes.len() < expected {
            return Err(format!(
                "File is {got} bytes but expected at least {expected}",
                got = bytes.len()
            ));
        }

        // ── Decode arrays ───────────────────────────────────────────────────
        let log_mass = read_f64_slice(&bytes[4..4 + 8 * n]);
        let conc_raw = read_f64_slice(&bytes[4 + 8 * n..4 + 16 * n]);

        // Store concentrations in log space for log-log interpolation.
        let log_conc: Vec<f64> = conc_raw.iter().map(|c| c.log10()).collect();

        Ok(Self { log_mass, log_conc })
    }

    /// Query the concentration at a given halo mass.
    ///
    /// * `mass_msun_h` – halo mass in Msun/h (linear, *not* log).
    ///
    /// Returns `None` if the mass falls outside the tabulated range.
    pub fn concentration(&self, mass_msun_h: f64) -> Option<f64> {
        let lm = mass_msun_h.log10();
        let log_c = self.interp_log_conc(lm)?;
        Some(10_f64.powf(log_c))
    }

    /// Same as [`concentration`] but accepts log10(M / [Msun/h]) directly.
    pub fn concentration_from_log_mass(&self, log_mass: f64) -> Option<f64> {
        let log_c = self.interp_log_conc(log_mass)?;
        Some(10_f64.powf(log_c))
    }

    /// Returns the tabulated log10(M) range as `(min, max)`.
    pub fn log_mass_range(&self) -> (f64, f64) {
        (
            *self.log_mass.first().unwrap(),
            *self.log_mass.last().unwrap(),
        )
    }

    // ── Internal ─────────────────────────────────────────────────────────────

    /// Linear interpolation in log10(M) – log10(c) space.
    fn interp_log_conc(&self, lm: f64) -> Option<f64> {
        let n = self.log_mass.len();
        let lo = *self.log_mass.first().unwrap();
        let hi = *self.log_mass.last().unwrap();

        if lm < lo || lm > hi {
            return None;
        }

        // Binary search for the bracketing interval.
        let idx = self
            .log_mass
            .partition_point(|&x| x <= lm)
            .saturating_sub(1)
            .min(n - 2);

        let x0 = self.log_mass[idx];
        let x1 = self.log_mass[idx + 1];
        let y0 = self.log_conc[idx];
        let y1 = self.log_conc[idx + 1];

        let t = (lm - x0) / (x1 - x0);
        Some(y0 + t * (y1 - y0))
    }
}

/// Decode a byte slice into a `Vec<f64>` assuming little-endian f64 values.
fn read_f64_slice(bytes: &[u8]) -> Vec<f64> {
    bytes
        .chunks_exact(8)
        .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
        .collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;

    /// Smoke-test: load the table and check a midpoint query is in a
    /// physically plausible range.  Run with:
    ///   cargo test -- --nocapture
    #[test]
    fn test_load_and_query() {
        let table = MassConcentrationTable::from_file("concentration_table.bin")
            .expect("Failed to load table – did you run gen_concentration_table.py?");

        let (lm_min, lm_max) = table.log_mass_range();
        println!("log10(M) range: [{lm_min}, {lm_max}]");

        // Milky-Way-ish halo at ~10^12 Msun/h
        let c = table.concentration(1e12).expect("Mass in range");
        println!("c(10^12 Msun/h) = {c:.3}");
        assert!(
            (3.0..30.0).contains(&c),
            "Concentration {c} outside plausible range"
        );

        // Out-of-range should return None
        assert!(table.concentration(1e5).is_none());
        assert!(table.concentration(1e15).is_none());
    }
}
