use std::{f64::consts::PI, path::MAIN_SEPARATOR_STR, str::FromStr};

use crate::hydrostatics::GG;

const M_IN_KPC: f64 = 3.0857e19;
const KG_IN_MSUN: f64 = 1.989e30;

const HH: f64 = 0.671; // normalization constant
const HUBBLE: f64 = 0.0671; // hubble const km / s kpc
const _RHOCRIT: f64 = 1.8791 * (HH * HH) * 1e-26 * (M_IN_KPC * M_IN_KPC * M_IN_KPC) / KG_IN_MSUN; // Halobos calculation, contains magic numbers I don't understand
const RHOCRIT: f64 = 3.0 * (HUBBLE * HUBBLE) / (8.0 * PI * GG); // My calculation

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

    pub fn r200(&self) -> Result<f64, String> {
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

    pub fn m200(&self) -> Result<f64, String> {
        Ok(self.enclosed_mass(self.r200()?))
    }

    pub fn c200(&self) -> Result<f64, String> {
        match self {
            Halo::NFW(_rho_s, r_s) => Ok(self.r200()? / r_s),
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
