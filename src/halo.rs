use std::path::Path;
use std::sync::OnceLock;
use std::{f64::consts::PI, str::FromStr};

use crate::concentration_table::*;
use crate::constants::*;

pub enum Halo {
    NFW(f64, f64), //rho_s, r_s
}

impl Halo {
    pub fn scale_radius(&self) -> f64 {
        match self {
            Halo::NFW(_rho_s, r_s) => *r_s,
        }
    }

    pub fn enclosed_mass(&self, r: f64) -> f64 {
        if r < 0.0 {
            return 0.0;
        }
        match self {
            Halo::NFW(rho_s, r_s) => {
                4.0 * PI * rho_s * r_s.powi(3) * (((r_s + r) / r_s).ln() - (r / (r_s + r)))
            }
        }
    }

    pub fn r_crit(&self) -> f64 {
        match self {
            Halo::NFW(rho_s, r_s) => {
                // for c = rho_s / rho_crit
                // solving r*(r_s + r)^2 = r_s^3 c
                let c = rho_s / RHOCRIT;
                let factor = (3.0
                    * 3.0_f64.sqrt()
                    * (27.0 * r_s.powi(6) * c.powi(2) + 4.0 * r_s.powi(6) * c).sqrt()
                    + 27.0 * r_s.powi(3) * c
                    + 2.0 * r_s.powi(3))
                .cbrt();
                ((factor / 2.0_f64.cbrt()) + (2.0_f64.cbrt() * r_s.powi(2) / factor) - 2.0 * r_s)
                    / 3.0
            }
        }
    }

    pub fn r200(&self) -> Result<f64, String> {
        match self {
            Halo::NFW(_rho_s, r_s) => {
                let c200 = self.c200()?;
                Ok(c200 * r_s)
            }
            #[allow(unreachable_patterns)]
            _ => {
                dbg!("hehe illegal code is running");
                let mut r200 = self.scale_radius(); // initial guess (definitely wrong but right order of magnitude)
                let mut rho_enc = self.enclosed_mass(r200) / (4.0 * PI * r200.powi(3) / 3.0);
                let mut diff = (rho_enc - 200.0 * RHOCRIT) / (200.0 * RHOCRIT);

                let mut n = 0;
                const TOLERANCE: f64 = 1e-5;
                while diff.abs() > TOLERANCE {
                    r200 += diff * self.scale_radius() * (1.0_f64 - 1e-3).powi(n);
                    rho_enc = self.enclosed_mass(r200) / (4.0 * PI * r200.powi(3) / 3.0);
                    diff = (rho_enc - 200.0 * RHOCRIT) / (200.0 * RHOCRIT);
                    n += 1;

                    if n >= 1000 {
                        return Err(String::from_str("Failed to converge on r200!").unwrap());
                    }
                }
                Ok(r200)
            }
        }
    }

    pub fn m200(&self) -> Result<f64, String> {
        Ok(self.enclosed_mass(self.r200()?))
    }

    pub fn c200(&self) -> Result<f64, String> {
        match self {
            Halo::NFW(rho_s, _r_s) => {
                let mut c200: f64 = 1.0; // initial guess, def wrong, right order of magnitude
                let mut value = ((1.0 + c200).ln() - (c200 / (1.0 + c200))) / c200.powi(3);
                let target = (1.0 / 3.0) * 200.0 * RHOCRIT / rho_s;
                let mut diff = (value - target) / target;

                let mut n = 0;
                const TOLERANCE: f64 = 1e-5;
                while diff.abs() > TOLERANCE {
                    c200 += diff * (1.0_f64 - 1e-3).powi(n);
                    value = ((1.0 + c200).ln() - (c200 / (1.0 + c200))) / c200.powi(3);
                    diff = (value - target) / target;
                    n += 1;

                    if n >= 1000 {
                        return Err(format!(
                            "Failed to converge on c200!\ndiff = {}, target = {}, value = {}",
                            diff, target, value
                        ));
                    }
                }
                Ok(c200)
            }
        }
    }

    pub fn deviation(&self) -> Result<f64, String> {
        match self {
            Halo::NFW(_rho_s, _r_s) => {
                // https://arxiv.org/abs/1402.7073 eqn 8
                let sigma = 0.11;
                let c200 = self.c200()?;
                let m200 = self.m200()?;
                Ok((c200.log10() + 0.101 * (m200 * HH / 1e12).log10() - 0.905) / sigma)
            }
        }
    }
}

pub fn m200_c200_to_rs_rhos(m200: f64, c200: f64) -> (f64, f64) {
    let mass_integral_factor = (1.0 + c200).ln() - (c200 / (1.0 + c200));
    let rho_s = 200.0 * RHOCRIT * (c200.powi(3) / (3.0 * mass_integral_factor));
    let geometric_factor = 4.0 * PI * c200.powi(3) / 3.0;
    let r_s = (m200 / (geometric_factor * 200.0 * RHOCRIT)).cbrt();

    (r_s, rho_s)
}

pub fn rs_rhos_to_m200_c200(r_s: f64, rho_s: f64) -> (f64, f64) {
    let halo = Halo::NFW(rho_s, r_s);
    let c200 = halo.c200().unwrap();
    let m200 = halo.enclosed_mass(r_s * c200);

    (m200, c200)
}

#[derive(Clone, Debug)]
pub enum McrSource {
    DuttonMaccio2014,
    DiemerJoyce2019,
}

static DIEMER_JOYCE_TABLE: OnceLock<MassConcentrationTable> = OnceLock::new();

pub fn init_diemer_joyce<P: AsRef<Path>>(path: P) -> Result<(), String> {
    // OnceLock::set fails if already initialised — that's fine, just ignore it.
    let _ = DIEMER_JOYCE_TABLE.set(MassConcentrationTable::from_file(path)?);
    Ok(())
}

pub fn deviation(m200: f64, c200: f64, source: McrSource) -> f64 {
    let (median_log10c, sigma_log10c) = mass_concentration_relation(m200, source);
    (c200.log10() - median_log10c) / sigma_log10c
}

pub fn mass_concentration_relation(m200: f64, source: McrSource) -> (f64, f64) {
    match source {
        McrSource::DuttonMaccio2014 => {
            // Dutton & Maccio 2014, z=0

            let a = 0.905;
            let b = -0.101;
            let sigma_log10c = 0.11; // intrinsic scatter

            let log10m = (m200 * HH / 1e12).log10();
            let mean_log10c = a + b * log10m;

            (mean_log10c, sigma_log10c)
        }
        McrSource::DiemerJoyce2019 => {
            // m200 is in Msun; the table expects Msun/h.
            let mass_msun_h = m200 / HH;
            let c = diemer_joyce_table()
                .concentration(mass_msun_h)
                .unwrap_or_else(|| {
                    panic!(
                        "DiemerJoyce2019: mass {:.3e} Msun/h is outside the table range",
                        mass_msun_h
                    )
                });
            let median_log10c = c.log10();
            (median_log10c, DIEMER_JOYCE_SCATTER)
        }
    }
}

#[inline]
fn diemer_joyce_table() -> &'static MassConcentrationTable {
    DIEMER_JOYCE_TABLE
        .get()
        .expect("DiemerJoyce2019: call init_diemer_joyce(path) before querying")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() {
        // Idempotent — safe to call from multiple tests.
        init_diemer_joyce("concentration_table.bin").expect("load table");
    }

    #[test]
    fn rs_rhos_to_m200_c200_conversion() {
        let c200: f64 = rand::random_range(0.0..20.0);
        let m200: f64 = (rand::random_range(7.0..15.0) as f64).exp();
        dbg!((m200, c200));

        let (rs, rhos) = m200_c200_to_rs_rhos(m200, c200);

        let (new_m200, new_c200) = rs_rhos_to_m200_c200(rs, rhos);

        let m200_err = (new_m200 - m200) / m200;
        let c200_err = (new_c200 - c200) / c200;

        if m200_err > 1e-4 || c200_err > 1e-4 {
            panic!("m200 error: {m200_err:.5}, c200 error: {c200_err:.5}");
        }
    }

    #[test]
    fn test_r_crit() {
        let halo = Halo::NFW(240.0 * RHOCRIT, 4.2);
        let r_crit = halo.r_crit();
        assert!((r_crit - 23.3785).abs() < 1e-4);
    }

    #[test]
    fn diemer_joyce_median_plausible() {
        setup();
        let (log10c, sigma) = mass_concentration_relation(1e12, McrSource::DiemerJoyce2019);
        let c = 10_f64.powf(log10c);
        println!("c(10^12 Msun) = {c:.3},  sigma = {sigma}");
        assert!((3.0..20.0).contains(&c));
        assert_eq!(sigma, DIEMER_JOYCE_SCATTER);
    }

    #[test]
    fn deviation_at_median_is_zero() {
        setup();
        let m200 = 1e11; // Msun
        let (log10c, _) = mass_concentration_relation(m200, McrSource::DiemerJoyce2019);
        let c_median = 10_f64.powf(log10c);
        let dev = deviation(m200, c_median, McrSource::DiemerJoyce2019);
        assert!(
            dev.abs() < 1e-10,
            "deviation at median should be ~0, got {dev}"
        );
    }

    #[test]
    fn dutton_maccio_unchanged() {
        // Regression guard — make sure the DM14 arm is untouched.
        let (log10c, sigma) = mass_concentration_relation(1e12, McrSource::DuttonMaccio2014);
        let c = 10_f64.powf(log10c);
        println!("DM14 c(10^12 Msun) = {c:.3}");
        assert!((3.0..20.0).contains(&c));
        assert_eq!(sigma, 0.11);
    }
}
