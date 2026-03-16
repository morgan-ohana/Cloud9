use std::f64::consts::PI;

// System Values
pub const RHO_IGM: f64 = 133.02; // mean density of universe M_sun / kpc^3
pub const UVB_TEMP: f64 = 1e4; //FIND REAL VALUE

// MCR
pub const DIEMER_JOYCE_SCATTER: f64 = 0.16; // intrinsic scatter in log10(c), Diemer & Joyce 2019

// Physical Constants
pub const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
pub const M_PROTON: f64 = 8.41e-58; // Proton mass in Msun
pub const MP_OVER_KB: f64 = 1.15349467e35 * MOLECULAR_WEIGHT; // Particle mass over boltzmann constant s^2 K / kpc^2
pub const MOLECULAR_WEIGHT: f64 = 0.5;
pub const DISTANCE: f64 = 5e3; // 5 MPC or 5000 KPC
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

const _M_IN_KPC: f64 = 3.0857e19;
const _KG_IN_MSUN: f64 = 1.989e30;
const _RHOCRIT: f64 =
    1.8791 * (HH * HH) * 1e-26 * (_M_IN_KPC * _M_IN_KPC * _M_IN_KPC) / _KG_IN_MSUN; // Halobos calculation, contains magic numbers I don't understand
