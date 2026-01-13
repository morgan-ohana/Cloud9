use std::f64::consts::PI;

const SPACIAL_GRID_NUM: usize = 10000;
const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
const M_P: f64 = 8.41e-58

fn get_r_points(bounds:(f64, f64)) -> Vec<f64> {
    let mut r_points = Vec::with_capacity(SPACIAL_GRID_NUM);
    for i in 0..SPACIAL_GRID_NUM - 1 {
        r_points.push((bounds.0.ln() + (i as f64) * (bounds.1.ln() - bounds.0.ln()) / ((SPACIAL_GRID_NUM - 1) as f64)).exp())
    }
    r_points
}

fn get_rho_points<T: Fn(f64)->f64>(rho: T, r_points: Vec<f64>) -> Vec<f64> {
    let mut rho_points = Vec::with_capacity(r_points.len());
    for i in 0..r_points.len() - 1 {
        rho_points.push(rho(r_points[i]))
    }
    rho_points
}

fn get_force_points(rho_points: Vec<f64>, r_points: Vec<f64>) -> Vec<f64> {
    let mut force_points = Vec::with_capacity(r_points.len());
    let mut enclosed_mass = 0.0;
    for i in 1..r_points.len() - 1 {
        let vol = (4.0 * PI / 3.0) * (r_points[i].powi(3) - rho_points[i-1].powi(3));
        let ave_rho = (rho_points[i] + rho_points[i-1]) / 2.0;
        enclosed_mass += ave_rho * vol;
        force_points.push( - GG * enclosed_mass / r_points[i].powi(2));
    }
    force_points
}

fn get_hydrostatic_profile(r_points: Vec<f64>, external_field: Vec<f64>, temperature_points: Vec<f64>, rho_center: f64) -> Vec<f64> {
    let mut rho_points = Vec::with_capacity(r_points.len());
    rho_points.push(rho_center);

    // For hydrostatics f = 0 = dP + f_grav dr + f_ext dr
    // We also know for an ideal gas P = rho KT/m so we have:
    // drho = - (m/kT) * (f_grav + f_ext) dr

    let mut enclosed_mass = 0.0;
    for i in 1..r_points.len() - 1 {
        let dr = r_points[i] - r_points[i-1];
        
        // f_ext dr via trapezoid
        let external_piece = dr * (external_field[i] + external_field[i-1]) / 2.0
    
        // enclosed mass 
        let vol = (4.0 * PI / 3.0) * (r_points[i-1].powi(3) - rho_points[i-2].powi(3));
        enclosed_mass += vol * (rho_points[i-1] + rho_points[i-2]) / 2.0;
        let f_grav = - GG * enclosed_mass / r_points[i].powi(2);

        let thermo_prefactor = ((M_P/temperature_points[i]) + (M_P/temperature_points[i-1])) / 2.0;
        let drho = - thermo_prefactor * (external_piece + f_grav * dr);

        rho_points.push(rho_points.last().unwrap() + drho);
    }

    rho_points
}