use std::f64::consts::PI;

use crate::plotting::plot_function;

const SPACIAL_GRID_NUM: usize = 1000;
pub const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
pub const KM_IN_KPC: f64 = 3.086e16;
pub const CM_IN_KPC: f64 = 3.086e21;
//const K_B: f64 = 7.29e-93; // Boltzmanns constant Msun kpc^2 / s^2 K
pub const MP_OVER_KB: f64 = 1.15349467e35 * MOLECULAR_WEIGHT; // Particle mass over boltzmann constant s^2 K / kpc^2
pub const M_PROTON: f64 = 8.41e-58; // Proton mass in Msun
pub const MOLECULAR_WEIGHT: f64 = 0.5;
pub const DISTANCE:f64 = 5e3; // 5 MPC or 5000 KPC
pub const ARC_MIN: f64 = PI / 10800.0;

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

fn argswhere(vec: &Vec<f64>, value: f64) -> (usize, usize) {
    let mut low = 0;
    let mut high = vec.len() - 1;

    while high - low > 1 {
        let mid = (high + low) / 2;
        match value > vec[mid] {
            true => low = mid,
            false => high = mid,
        }
    }

    (low, high)
}

fn get_column_density(rho_points: Vec<f64>, r_points: &Vec<f64>) -> Vec<f64> {
    // Σ(R) = 2 ∫_0^∞ ρ(√(z² + R²)) dz = 2 ∫_R^∞ [ρ(r) * r / √(r² - R²)] dr

    let mut sigma = vec![0.0; r_points.len()];
    let r_max = *r_points.last().unwrap();
    let z_points = r_points.clone();

    for i in 0..r_points.len() {
        let r_projected = r_points[i];

        for j in 1..z_points.len() {
            let z1 = z_points[j - 1];
            let z2 = z_points[j];
            let dz = z2 - z1;
            let r1 = (z1.powi(2) + r_projected.powi(2)).sqrt();
            let r2 = (z2.powi(2) + r_projected.powi(2)).sqrt();

            if r1 > r_max {
                break;
            }

            let rho1 = {
                let (low, high) = argswhere(&r_points, r1);
                let t = (r1 - r_points[low]) / (r_points[high] - r_points[low]);
                t * rho_points[high] + (1.0 - t) * rho_points[low]
            };
            let rho2 = {
                let (low, high) = argswhere(&r_points, r2);
                let t = (r2 - r_points[low]) / (r_points[high] - r_points[low]);
                t * rho_points[high] + (1.0 - t) * rho_points[low]
            };
            sigma[i] += (rho1 + rho2) * dz;
        }
    }
    sigma
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
    bounds: (f64, f64),
    rho_center: f64,
) {
    let rho = |r: f64| -> f64 {
        rho_s / ((r / r_s).powf(gamma) * (1.0 + (r / r_s).powf(alpha)).powf((beta - gamma) / alpha))
    };

    let r_points = get_r_points(bounds);

    let dark_matter_rho_points = get_rho_points(rho, &r_points);

    let external_field = get_force_points(dark_matter_rho_points, &r_points);

    let temperature_points = vec![temp; r_points.len()];

    let rho_points =
        get_hydrostatic_profile(&r_points, external_field, temperature_points, rho_center);

    let number_density = {
        let mut num_density = Vec::with_capacity(rho_points.len());
        for i in 0..rho_points.len() {
            num_density.push(rho_points[i] / M_PROTON)
        }
        num_density
    };

    plot_function(
        &r_points,
        &number_density,
        "profile.png",
        "hydrostatic profile",
        "r (kpc)",
        "n_H (num / kpc^3)",
    )
    .unwrap();

    let column_density = {
        let mut col_dens = get_column_density(number_density, &r_points);
        for i in 0..col_dens.len() {
            col_dens[i] /= CM_IN_KPC.powi(2);
        }
        col_dens
    };

    let angular_points = {
        let mut ang_points = r_points.clone();
        for i in 0..ang_points.len() {
            ang_points[i] /= DISTANCE; // radians
            ang_points[i] /= ARC_MIN; // arc mins
        }
        ang_points
    };

    plot_function(
        &angular_points,
        &column_density,
        "column.png",
        "hydrostatic column density",
        "r (arcmin)",
        "n_H (num / cm^2)",
    )
    .unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_column_density_inv() {
        let r_points = get_r_points((1e-3, 1e3));
        let rho_points = get_rho_points(|r: f64| r.powi(-2), &r_points);
        //dbg!(&rho_points);
        let column_density = get_column_density(rho_points, &r_points);

        let r_max = r_points.last().unwrap();
        let analytic_column_density = get_rho_points(
            |r: f64| 2.0 * ((r_max.powi(2) - r.powi(2)).sqrt() / r).atan() * r.powi(-1),
            &r_points,
        );

        let err = {
            let mut err = Vec::with_capacity(r_points.len());
            for i in 0..r_points.len() {
                err.push(
                    (column_density[i] - analytic_column_density[i])
                        / (analytic_column_density[i] + 1e-15),
                )
            }
            err
        };

        let rms_err: f64 = {
            let mut rms_err: f64 = 0.0;
            for i in 0..err.len() - 1 {
                rms_err += err[i].powi(2) * (r_points[i+1] - r_points[i])
            }
            rms_err /= r_points.last().unwrap();
            rms_err.sqrt()
        };

        if rms_err > 0.01 || !rms_err.is_finite() {
            //dbg!(err);
            plot_function(&r_points, &err, "column_err_check.png", "Column Density Error", "r (kpc)", "% err").unwrap();
            panic!("rms error too high! rms_err = {rms_err}")
        }
    }
    
    #[test]
    fn test_column_density_lin() {
        let r_points = get_r_points((1e-3, 1e3));
        let rho_points = get_rho_points(|_r: f64| 1.0, &r_points);
        //dbg!(&rho_points);
        let column_density = get_column_density(rho_points, &r_points);

        let r_max = r_points.last().unwrap();
        let analytic_column_density = get_rho_points(
            |r: f64| 2.0 * (r_max.powi(2) - r.powi(2)).sqrt(),
            &r_points,
        );

        let err = {
            let mut err = Vec::with_capacity(r_points.len());
            for i in 0..r_points.len() {
                err.push(
                    (column_density[i] - analytic_column_density[i])
                        / (analytic_column_density[i] + 1e-15),
                )
            }
            err
        };

        let rms_err: f64 = {
            let mut rms_err: f64 = 0.0;
            for i in 0..err.len() {
                rms_err += err[i].powi(2)
            }
            rms_err /= err.len() as f64;
            rms_err.sqrt()
        };

        if rms_err > 0.01 || !rms_err.is_finite() {
            //dbg!(err);
            plot_function(&r_points, &err, "column_err_check2.png", "Column Density Error", "r (kpc)", "% err").unwrap();
            panic!("rms error too high! rms_err = {rms_err}")
        }
    }

}
