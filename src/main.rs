use std::f64::consts::PI;

use crate::hydrostatics::isothermal_abg_background;

mod hydrostatics;
mod plotting;

fn main() {
    let target_init_num_col_density = 1e19;
    let target_init_col_density: f64 = target_init_num_col_density * hydrostatics::M_PROTON * hydrostatics::CM_IN_KPC.powi(2);
    let temperature = 1e4;
    let mut sound_speed_squared = temperature / (hydrostatics::MP_OVER_KB); // kpc^2 / s^2
    sound_speed_squared *= hydrostatics::KM_IN_KPC.powi(2); // km^2 / s^2
    
    // r_s = c_s / sqrt(4pi G rho_c)
    // uniform sphere approx:
    // Sigma ~ 2*r_s*rho_c = 2 c_s sqrt(rho_c / 4pi G) => rho_c ~ 4pi G * (Sigma / 2c_s)^2 = pi G Sigma^2 / c_s^2 
    // rho_c (kpc^-3) ~ pi G * (n_sigma (cm^-2) * MP * CMinKPC^2)^2 / c_s^2 
    // Units: [rho_c] = M_sun kpc^-3 = [G * sigma^2 / c_s^2] = km^2 kpc M_sun^-1 s^-2 * M_sun^2 kpc^-4 * [c_s]^-2 = km^2 kpc^-3 M_sun s^-2 * [c_s]^-2
    // => [c_s^2] = km^2 kpc^-3 M_sun s^-2 / M_sun kpc^-3 = km^2 / s^2
    
    let rho_center_approx = PI * hydrostatics::GG * (target_init_col_density).powi(2) / sound_speed_squared;
    let scale_radius = sound_speed_squared.sqrt() / (4.0 * PI * hydrostatics::GG * rho_center_approx);
    println!("Center rho set to: {rho_center_approx}");

    isothermal_abg_background(
        temperature,
        1.0,
        3.0,
        1.0,
        0.0e7,
        3.0,
        (1e-1, 10.0*scale_radius),
        rho_center_approx,
    );
}
