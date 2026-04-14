use std::f64::consts::PI;

// System Values
pub const INNER_BOUND: f64 = 0.015;
pub const RHO_IGM: f64 = 133.02; // mean density of universe M_sun / kpc^3
pub const DISTANCE: f64 = 4.6e3; // 4.6 MPC or 4660 KPC https://iopscience.iop.org/article/10.3847/1538-4357/acdcf5/pdf

// MCR
pub const DIEMER_JOYCE_SCATTER: f64 = 0.16; // intrinsic scatter in log10(c), Diemer & Joyce 2019

// Physical Constants
pub const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
pub const M_PROTON: f64 = 8.41e-58; // Proton mass in Msun
pub const MP_OVER_KB: f64 = 1.15349467e35; // Proton mass over boltzmann constant s^2 K / kpc^2
pub const HH: f64 = 0.671; // normalization factor for hubble constant
pub const HUBBLE: f64 = 0.0671; // hubble const km / s kpc
pub const RHOCRIT: f64 = 3.0 * (HUBBLE * HUBBLE) / (8.0 * PI * GG); // My calculation (124.927 M_sun / kpc^3)

//Conversions
pub const KM_IN_KPC: f64 = 3.086e16;
pub const CM_IN_KPC: f64 = 3.086e21;
//const K_B: f64 = 7.29e-93; // Boltzmanns constant Msun kpc^2 / s^2 K
pub const ARC_MIN: f64 = PI / 10800.0;
pub const S_IN_GYR: f64 = 3.154e16;
pub const G_IN_MSUN: f64 = 1.988e33;
