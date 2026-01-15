use std::f64::consts::PI;

use crate::plotting::plot_function;

const SPACIAL_GRID_NUM: usize = 100000;
const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
const KM_IN_KPC: f64 = 3.086e16;
const K_B: f64 = 7.29e-93; // Boltzmanns constant Msun kpc^2 / s^2 K
const MP_OVER_KB: f64 = 1.15349467e35; // Proton mass over boltzmann constant s^2 K / kpc^2

fn get_r_points(bounds: (f64, f64)) -> Vec<f64> {
    let mut r_points = Vec::with_capacity(SPACIAL_GRID_NUM);
    for i in 0..SPACIAL_GRID_NUM {
        r_points.push(
            (bounds.0.ln()
                + (i as f64) * (bounds.1.ln() - bounds.0.ln()) / ((SPACIAL_GRID_NUM - 1) as f64))
                .exp(),
        )
    }
    r_points
}

fn get_rho_points<T: Fn(f64) -> f64>(rho: T, r_points: &Vec<f64>) -> Vec<f64> {
    let mut rho_points = Vec::with_capacity(r_points.len());
    for i in 0..r_points.len() {
        rho_points.push(rho(r_points[i]));
        //dbg!(i);
        //dbg!(rho_points.len());
    }
    rho_points
}

fn get_force_points(rho_points: Vec<f64>, r_points: &Vec<f64>) -> Vec<f64> {
    let mut force_points = Vec::with_capacity(r_points.len());
    force_points.push(0.0);
    let mut enclosed_mass = 0.0;
    for i in 1..r_points.len() {
        let vol = (4.0 * PI / 3.0) * (r_points[i].powi(3) - r_points[i - 1].powi(3));
        let ave_rho = (rho_points[i] + rho_points[i - 1]) / 2.0;
        enclosed_mass += ave_rho * vol;
        force_points.push(-GG * enclosed_mass / r_points[i].powi(2)); // [km^2 kpc / M_sun s^2] * Msun / kpc^2 = [km^2 / s^2] / kpc
        force_points[i] /= KM_IN_KPC * KM_IN_KPC; // kpc / s^2 
    }
    force_points
}

fn get_hydrostatic_profile(
    r_points: &Vec<f64>,
    external_field: Vec<f64>,
    temperature_points: Vec<f64>,
    rho_center: f64,
) -> Vec<f64> {
    let mut rho_points = Vec::with_capacity(r_points.len());
    rho_points.push(rho_center);
    rho_points.push(rho_center);

    // For hydrostatics f = 0 = -dP + f_grav dr + f_ext dr
    // We also know for an ideal gas P = rho KT/m so we have:
    // drho = (m/kT) * rho * (a_grav + a_ext) dr

    // Units: rho in M_sun / kpc^3 => [P] = [rho KT/m] = M_sun/ kpc^3 * kpc^2 / s^2 = (M_sun kpc / s^2) / kpc^2

    let mut enclosed_mass = 0.0;
    for i in 2..r_points.len() {
        let dr = r_points[i] - r_points[i - 1];

        // f_ext dr via trapezoid
        let external_piece = dr * (external_field[i] + external_field[i - 1]) / 2.0;

        // enclosed mass
        let vol = (4.0 * PI / 3.0) * (r_points[i - 1].powi(3) - r_points[i - 2].powi(3));
        enclosed_mass += vol * (rho_points[i - 1] + rho_points[i - 2]) / 2.0;
        let mut f_grav = -GG * enclosed_mass / r_points[i].powi(2); // [km^2 kpc / M_sun s^2] * Msun / kpc^2 = [km^2 / s^2] / kpc
        f_grav /= KM_IN_KPC * KM_IN_KPC; // kpc / s^2

        let thermo_prefactor =
            ((MP_OVER_KB / temperature_points[i]) + (MP_OVER_KB / temperature_points[i - 1])) / 2.0;
        let drho = thermo_prefactor * rho_points[i - 1] * (external_piece + f_grav * dr);
        //dbg!(drho);

        rho_points.push((rho_points.last().unwrap() + drho).max(0.0));
    }

    rho_points
}

pub fn isothermal_abg_background(
    temp: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    rho_s: f64,
    r_s: f64,
    rho_center: f64,
) {
    let rho = |r: f64| -> f64 {
        rho_s / ((r / r_s).powf(gamma) * (1.0 + (r / r_s).powf(alpha)).powf((beta - gamma) / alpha))
    };

    let r_points = get_r_points((r_s * 1e-3, r_s * 1e2));

    let dark_matter_rho_points = get_rho_points(rho, &r_points);

    let external_field = get_force_points(dark_matter_rho_points, &r_points);

    let temperature_points = vec![temp; r_points.len()];

    let rho_points =
        get_hydrostatic_profile(&r_points, external_field, temperature_points, rho_center);

    plot_function(
        &r_points,
        &rho_points,
        "profile.png",
        "hydrostatic profile",
        "r (kpc)",
        "rho (M_sun / kpc^3)",
    )
    .unwrap();
}
