use std::f64::consts::PI;

use crate::constants::*;
use crate::halo::Halo;
use crate::plotting::plot_functions;

const SPACIAL_GRID_NUM: usize = 1000;

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

fn get_total_mass(rho_points: &Vec<f64>, r_points: &Vec<f64>) -> f64 {
    let mut enclosed_mass = 0.0;
    for i in 1..r_points.len() {
        let vol = (4.0 * PI / 3.0) * (r_points[i].powi(3) - r_points[i - 1].powi(3));
        let ave_rho = (rho_points[i] + rho_points[i - 1]) / 2.0;
        enclosed_mass += ave_rho * vol;
    }
    enclosed_mass
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

fn get_hydrostatic_profile_outwards(
    r_points: &Vec<f64>,

    external_field: Vec<f64>,

    temperature_points: &Vec<f64>,

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

fn get_hydrostatic_profile_inwards(
    r_points: &Vec<f64>,
    dm_rho_points: Vec<f64>,
    temperature_points: Vec<f64>,
) -> Vec<f64> {
    let mut gas_rho_points = vec![0.0; SPACIAL_GRID_NUM];
    gas_rho_points[SPACIAL_GRID_NUM - 1] = RHO_IGM;
    gas_rho_points[SPACIAL_GRID_NUM - 2] = RHO_IGM;
    gas_rho_points[SPACIAL_GRID_NUM - 3] = RHO_IGM;

    // For hydrostatics f = 0 = -dP + f_grav dr + f_ext dr
    // We also know for an ideal gas P = rho KT/m so we have:
    // drho = (m/kT) * rho * (a_grav + a_ext) dr
    //
    // dP/dr = - G M_enc rho_g / r^2
    // (r^2/rho_g) dP/dr = - G M_enc
    // M_enc = 4pi ∫_0^r (rho) r'^2 dr' where rho = rho_dm + rho_gas
    // Differentiating we get the second order ODE:
    // (2r/rho_g) dP/dr + (r^2/rho_g) d^2P/dr^2 - (r^2/rho_g^2)(dP/dr)(drho_g/dr) = -4piG rho r^2
    // => d^2P/dr^2 = (dP/dr)(drho_g/dr)/rho_g - 2(dP/dr)/r - 4piG rho rho_g
    // remembering P = rho_gas kT/m
    // d^2rho_g/dr^2 = (drho_g/dr)^2/rho_g - 2(drho_g/dr)/r - (4piGm/kT)rho * rho_g

    // Units: rho in M_sun / kpc^3 => [P] = [rho KT/m] = M_sun/ kpc^3 * kpc^2 / s^2 = (M_sun kpc / s^2) / kpc^2
    // [d^2rho_g/dr^2] = M_sun / kpc^5 and [rho^2] = M_sun^2 / kpc^6
    // => [Gm/kT] = kpc / M_sun
    // [G] = km^2 kpc / M_sun s^2 and [m/k] = s^2 K / kpc^2 and [T] = K
    // => [Gm/kT] = km^2 / M_sun kpc = kpc / M_sun (km/kpc)^2

    dbg!(r_points);
    let mut dgas_rho_dr: f64 = 0.0;
    for i in (1..SPACIAL_GRID_NUM - 3).rev() {
        dbg!(i);
        dbg!(r_points[i]);
        let force_prefactor = {
            let mut prefactor = 4.0 * PI * GG * MP_OVER_KB / temperature_points[i];
            prefactor /= KM_IN_KPC * KM_IN_KPC; // kpc / M_sun
            prefactor
        };

        let dr = r_points[i + 1] - r_points[i];

        let d2gas_rho_dr2 = {
            dbg!(dm_rho_points[i]);
            let rho_g = gas_rho_points[i + 1];
            let r = r_points[i];
            let term_1 = dgas_rho_dr.powi(2) / rho_g;
            dbg!(term_1);
            let term_2 = -2.0 * dgas_rho_dr / r;
            dbg!(term_2);
            let rho = rho_g + dm_rho_points[i]; //gas_rho_points[i] does not exist yet, use one point higher
            let term_3 = -1.0 * force_prefactor * rho * rho_g;
            dbg!(term_3);
            term_1 + term_2 + term_3
        };
        dbg!(d2gas_rho_dr2);

        dgas_rho_dr -= d2gas_rho_dr2 * dr; // -= because we are integrating inwards
        dbg!(dgas_rho_dr);

        gas_rho_points[i] = gas_rho_points[i + 1] - dgas_rho_dr * dr; // - delta rho because integrating inwards
        dbg!(gas_rho_points[i]);

        if gas_rho_points[i] < 0.0 || !gas_rho_points[i].is_finite() {
            panic!("FUCK");
        }
    }
    // idx 0 never reached in intergal, just flatten to remove rho(0)=0 point
    gas_rho_points[0] = gas_rho_points[1];

    gas_rho_points
}

pub fn isothermal_abg_background(
    temp: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    rho_s: f64,
    r_s: f64,
    rho_c: Option<f64>,
    bounds: (f64, f64),
) {
    let rho = |r: f64| -> f64 {
        rho_s / ((r / r_s).powf(gamma) * (1.0 + (r / r_s).powf(alpha)).powf((beta - gamma) / alpha))
    };

    let r_points = get_r_points(bounds);

    let dark_matter_rho_points = get_rho_points(rho, &r_points);

    let temperature_points = vec![temp; r_points.len()];

    let rho_points = match rho_c {
        None => {
            get_hydrostatic_profile_inwards(&r_points, dark_matter_rho_points, temperature_points)
        }
        Some(rho_c) => {
            let external_field = get_force_points(dark_matter_rho_points, &r_points);
            get_hydrostatic_profile_outwards(&r_points, external_field, &temperature_points, rho_c)
        }
    };

    let number_density = {
        let mut num_density = Vec::with_capacity(rho_points.len());
        for i in 0..rho_points.len() {
            num_density.push(rho_points[i] / M_PROTON)
        }
        num_density
    };

    plot_functions(
        &r_points,
        &vec![number_density.clone()],
        "profile.png",
        "hydrostatic profile",
        "r (kpc)",
        "n_H (num / kpc^3)",
        vec![None],
        vec![false],
        None,
        None,
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

    plot_functions(
        &angular_points,
        &vec![column_density],
        "column.png",
        "hydrostatic column density",
        "r (arcmin)",
        "n_H (num / cm^2)",
        vec![None],
        vec![false],
        None,
        None,
    )
    .unwrap();
}

pub fn isothermal_core_collapse_background(
    temp: f64,
    rho_s_0: f64,
    r_s_0: f64,
    collapse_progress: f64,
    rho_c: Option<f64>,
    bounds: (f64, f64),
    plot: bool,
) -> (Vec<f64>, Vec<f64>) {
    let rho = parametic_core_collapse(r_s_0, rho_s_0, collapse_progress);

    let r_points = get_r_points(bounds);

    let dark_matter_rho_points = get_rho_points(rho, &r_points);

    let temperature_points = vec![temp; r_points.len()];

    let rho_points = match rho_c {
        None => {
            get_hydrostatic_profile_inwards(&r_points, dark_matter_rho_points, temperature_points)
        }
        Some(rho_c) => {
            let external_field = get_force_points(dark_matter_rho_points, &r_points);
            get_hydrostatic_profile_outwards(&r_points, external_field, &temperature_points, rho_c)
        }
    };

    let number_density = {
        let mut num_density = Vec::with_capacity(rho_points.len());
        for i in 0..rho_points.len() {
            num_density.push(rho_points[i] / M_PROTON)
        }
        num_density
    };

    /*
    let halo = Halo::NFW(rho_s, r_s);

    let r_crit = halo.r_crit();
    let mut high = r_points.len() - 1;
    let mut low = 0;
    while high - low > 1 {
        let mid = (high + low) / 2;

        if r_points[mid] > r_crit {
            high = mid;
        } else {
            low = mid;
        }
    }
    dbg!(low);

    println!("rho(r_crit): {}", rho_points[low]);
    */

    if plot {
        plot_functions(
            &r_points,
            &vec![number_density.clone()],
            "profile.png",
            "hydrostatic profile",
            "r (kpc)",
            "n_H (num / kpc^3)",
            vec![None],
            vec![false],
            None,
            None,
        )
        .unwrap();
    }

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

    if plot {
        plot_functions(
            &angular_points,
            &vec![column_density.clone()],
            "column.png",
            "hydrostatic column density",
            "r (arcmin)",
            "n_H (num / cm^2)",
            vec![None],
            vec![false],
            None,
            None,
        )
        .unwrap();
    }

    (angular_points, column_density)
}

fn parametic_core_collapse(r_s_0: f64, rho_s_0: f64, tau: f64) -> impl Fn(f64) -> f64 {
    // https://arxiv.org/pdf/2406.10753 eqn 1 & 2
    let rho_s = rho_s_0
        * (2.033 + 0.7381 * tau + 7.264 * tau.powi(5) - 12.73 * tau.powi(7)
            + 9.915 * tau.powi(9)
            + (1.0 - 2.033) * (tau + 0.001).ln() / (0.001_f64).ln());
    let r_s = r_s_0
        * (0.7178 - 0.1026 * tau + 0.2474 * tau.powi(2) - 0.4079 * tau.powi(3)
            + (1.0 - 0.7178) * (tau + 0.001).ln() / (0.001_f64).ln());
    let r_c = r_s_0
        * (2.555 * tau.sqrt() - 3.632 * tau + 2.131 * tau.powi(2) - 1.415 * tau.powi(3)
            + 0.4683 * tau.powi(4));

    move |r: f64| -> f64 {
        rho_s / (((r.powi(4) + r_c.powi(4)).sqrt().sqrt() / r_s) * (1.0 + (r / r_s)).powi(2))
    }
}

pub fn evolution_profile(
    anchor_params: &[f64; 4],
    data: &Vec<(f64, f64)>,
    data_y_err: &Vec<(f64, f64)>,
) {
    let halo = Halo::NFW(anchor_params[0], anchor_params[1]);
    let bounds = (0.01 * anchor_params[1], halo.r_crit());
    let r_points = get_r_points(bounds);
    let temperature_points = vec![UVB_TEMP; r_points.len()];

    fn get_profile_and_mass(
        params: &[f64; 4],
        r_points: &Vec<f64>,
        temperature_points: &Vec<f64>,
    ) -> (Vec<f64>, f64) {
        let dm_rho = parametic_core_collapse(params[1], params[0], params[2]);
        let dm_rho_points = get_rho_points(dm_rho, r_points);
        let external_field = get_force_points(dm_rho_points, r_points);
        let rho_gas_points = get_hydrostatic_profile_outwards(
            r_points,
            external_field,
            temperature_points,
            params[3],
        );
        let total_mass = get_total_mass(&rho_gas_points, r_points);

        let number_density = {
            let mut num_density = Vec::with_capacity(rho_gas_points.len());
            for i in 0..rho_gas_points.len() {
                num_density.push(rho_gas_points[i] / M_PROTON)
            }
            num_density
        };

        let column_density = {
            let mut col_dens = get_column_density(number_density, &r_points);
            for i in 0..col_dens.len() {
                col_dens[i] /= CM_IN_KPC.powi(2);
            }
            col_dens
        };

        (column_density, total_mass)
    }

    // get anchor fit
    let (anchor_gas_rho_points, total_mass) =
        get_profile_and_mass(anchor_params, &r_points, &temperature_points);

    // get evolution snapshots

    let tau_snapshots = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0];
    let mut gas_rho_points_vec = Vec::with_capacity(tau_snapshots.len() + 1);
    let mut legends = Vec::with_capacity(tau_snapshots.len() + 1);
    gas_rho_points_vec.push(anchor_gas_rho_points);
    legends.push(Some(format!(
        "Modern Cloud-9, tau = {:.2}",
        anchor_params[2]
    )));
    let mut snapshot_params = anchor_params.clone();
    for tau in tau_snapshots {
        snapshot_params[2] = tau;
        // Note I don't reset snapshot rho_c between snapshots as snapshots are likely to be closer to each other than the anchor point
        let (mut snapshot_gas_rho_points, mut snapshot_total_mass) =
            get_profile_and_mass(&snapshot_params, &r_points, &temperature_points);

        let mut percent_diff = (snapshot_total_mass - total_mass) / total_mass;

        while percent_diff.abs() >= 0.01 {
            dbg!(percent_diff);
            dbg!(snapshot_params[3]);
            snapshot_params[3] -= percent_diff * anchor_params[3];

            (snapshot_gas_rho_points, snapshot_total_mass) =
                get_profile_and_mass(&snapshot_params, &r_points, &temperature_points);

            percent_diff = (snapshot_total_mass - total_mass) / total_mass
        }

        gas_rho_points_vec.push(snapshot_gas_rho_points);
        legends.push(Some(format!("tau = {}", tau)));
    }

    let angular_points = {
        let mut ang_points = r_points.clone();
        for i in 0..ang_points.len() {
            ang_points[i] /= DISTANCE; // radians
            ang_points[i] /= ARC_MIN; // arc mins
        }
        ang_points
    };

    let mut dashed = vec![true; gas_rho_points_vec.len()];
    dashed[0] = false;

    plot_functions(
        &angular_points,
        &gas_rho_points_vec,
        "figures/evolution.svg",
        "cloud evol",
        "r (kpc)",
        "rho (M_sun kpc^-3)",
        legends,
        dashed,
        Some(data),
        Some(data_y_err),
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
                rms_err += err[i].powi(2) * (r_points[i + 1] - r_points[i])
            }
            rms_err /= r_points.last().unwrap();
            rms_err.sqrt()
        };

        if rms_err > 0.01 || !rms_err.is_finite() {
            //dbg!(err);
            plot_functions(
                &r_points,
                &vec![err],
                "column_err_check.png",
                "Column Density Error",
                "r (kpc)",
                "% err",
                vec![None],
                vec![false],
                None,
            )
            .unwrap();
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
        let analytic_column_density =
            get_rho_points(|r: f64| 2.0 * (r_max.powi(2) - r.powi(2)).sqrt(), &r_points);

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
            plot_functions(
                &r_points,
                &vec![err],
                "column_err_check2.png",
                "Column Density Error",
                "r (kpc)",
                "% err",
                vec![None],
                vec![false],
                None,
            )
            .unwrap();
            panic!("rms error too high! rms_err = {rms_err}")
        }
    }
}
