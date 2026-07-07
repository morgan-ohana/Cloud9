use std::f64::consts::PI;
use std::sync::OnceLock;

use plotters::style::{FontDesc, IntoFont};
use rayon::prelude::*;

use crate::halo::{Halo, rs_rhos_to_m200_c200};
use crate::plotting::plot_functions;
use crate::temperature::{self, TnRelation};
use crate::{constants::*, fitting};

const SPACIAL_GRID_NUM: usize = 2000;

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

fn get_spacial_integral(rho_points: &Vec<f64>, r_points: &Vec<f64>) -> f64 {
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

fn get_potential_points(rho_points: Vec<f64>, r_points: &Vec<f64>) -> Vec<f64> {
    let mut potential_points = Vec::with_capacity(r_points.len());
    potential_points.push(0.0);
    let mut enclosed_mass = 0.0;
    for i in 1..r_points.len() {
        let vol = (4.0 * PI / 3.0) * (r_points[i].powi(3) - r_points[i - 1].powi(3));
        let ave_rho = (rho_points[i] + rho_points[i - 1]) / 2.0;
        enclosed_mass += ave_rho * vol;
        potential_points.push(-GG * enclosed_mass / r_points[i]); // [km^2 kpc / M_sun s^2] * Msun / kpc = [km^2 / s^2]
        potential_points[i] /= KM_IN_KPC * KM_IN_KPC; // kpc^2 / s^2 
    }
    potential_points
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

fn get_column_density_at_points(
    rho_points: Vec<f64>,
    r_points: &Vec<f64>,
    projected_points: &Vec<f64>,
) -> Vec<f64> {
    // Σ(R) = 2 ∫_0^∞ ρ(√(z² + R²)) dz = 2 ∫_R^∞ [ρ(r) * r / √(r² - R²)] dr

    let mut sigma = vec![0.0; projected_points.len()];
    let r_max = *r_points.last().unwrap();
    let z_points = r_points.clone();

    for i in 0..projected_points.len() {
        let r_projected = projected_points[i];

        for j in 1..z_points.len() {
            let z1 = z_points[j - 1];
            let z2 = z_points[j];
            let dz = z2 - z1;
            let r1 = (z1.powi(2) + r_projected.powi(2)).sqrt();
            let r2 = (z2.powi(2) + r_projected.powi(2)).sqrt();

            if r2 > r_max {
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

            if r2 > r_max {
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

fn get_hydrostatic_profile_outwards<T: Fn(f64, f64) -> f64, V: Fn(f64, f64) -> f64>(
    r_points: &Vec<f64>,
    external_field: Vec<f64>,
    temperature: T,
    neutral_fraction: V,
    rho_center: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut neutral_rho_points = Vec::with_capacity(r_points.len());
    let mut rho_points = Vec::with_capacity(r_points.len());

    let central_temp = temperature(0.0, rho_center);
    let neutral_rho_center = neutral_fraction(rho_center, central_temp) * rho_center;

    rho_points.push(rho_center);
    rho_points.push(rho_center);

    neutral_rho_points.push(neutral_rho_center);
    neutral_rho_points.push(neutral_rho_center);

    // For hydrostatics f = 0 = -dP + f_grav dr + f_ext dr
    // So dP = rho * (a_grav + a_ext) dr

    // We also know for an ideal gas P = rho KT/m so we have:

    // Units: rho in M_sun / kpc^3 => [P] = [rho KT/m] = M_sun/ kpc^3 * kpc^2 / s^2 = (M_sun kpc / s^2) / kpc^2

    let mut enclosed_mass = 0.0;

    for i in 2..r_points.len() {
        // current pressure & temp
        let temp = temperature(r_points[i - 1], rho_points[i - 1]);
        let neutral_frac = neutral_fraction(rho_points[i - 1], temp);
        let molecular_weight = 1.0 / (2.0 - neutral_frac);
        let mut press = rho_points[i - 1] * temp / (MP_OVER_KB * molecular_weight);

        // Calculate dP
        let dr = r_points[i] - r_points[i - 1];

        // f_ext dr via trapezoid

        let external_piece = dr * (external_field[i] + external_field[i - 1]) / 2.0;

        // enclosed mass

        let vol = (4.0 * PI / 3.0) * (r_points[i - 1].powi(3) - r_points[i - 2].powi(3));

        enclosed_mass += vol * (rho_points[i - 1] + rho_points[i - 2]) / 2.0;

        let mut f_grav = -GG * enclosed_mass / r_points[i].powi(2); // [km^2 kpc / M_sun s^2] * Msun / kpc^2 = [km^2 / s^2] / kpc

        f_grav /= KM_IN_KPC * KM_IN_KPC; // kpc / s^2

        let dpress = rho_points[i - 1] * (external_piece + f_grav * dr);

        press = (press + dpress).max(0.0);

        let rho_new = if press <= 0.0 {
            0.0
        } else {
            assert!(temp > 0.0);
            press * MP_OVER_KB * molecular_weight / temp
        };

        rho_points.push(rho_new);
        let temp = temperature(r_points[i], rho_points[i]);
        let neutral_frac = neutral_fraction(rho_points[i], temp);
        neutral_rho_points.push(rho_points[i] * neutral_frac);
    }

    (neutral_rho_points, rho_points)
}

fn get_hydrostatic_profile_inwards<T: Fn(f64, f64) -> f64>(
    r_points: &Vec<f64>,
    dm_rho_points: Vec<f64>,
    temperature: T,
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
            let mut prefactor =
                4.0 * PI * GG * MP_OVER_KB / temperature(r_points[i], gas_rho_points[i + 1]);
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

pub fn abg_background<T: Fn(f64, f64) -> f64, V: Fn(f64, f64) -> f64>(
    temperature: T,
    neutral_fraction: V,
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

    let neutral_rho_points = match rho_c {
        None => get_hydrostatic_profile_inwards(&r_points, dark_matter_rho_points, temperature),
        Some(rho_c) => {
            let external_field = get_force_points(dark_matter_rho_points, &r_points);
            get_hydrostatic_profile_outwards(
                &r_points,
                external_field,
                &temperature,
                &neutral_fraction,
                rho_c,
            )
            .0
        }
    };

    let number_density = {
        let mut num_density = Vec::with_capacity(neutral_rho_points.len());
        for i in 0..neutral_rho_points.len() {
            let num_dens = neutral_rho_points[i] / M_PROTON; // nH kpc^-3
            num_density.push(num_dens / CM_IN_KPC.powi(3))
        }
        num_density
    };

    plot_functions(
        &r_points,
        &vec![number_density.clone()],
        "profile.png",
        "hydrostatic profile",
        "r (kpc)",
        "n_H (num / cm^3)",
        vec![None],
        ("sans-serif", 12).into_font(),
        vec![false],
        None,
        None,
    )
    .unwrap();

    let column_density = {
        let mut col_dens = get_column_density(number_density, &r_points); // nH cm^-3 kpc
        for i in 0..col_dens.len() {
            col_dens[i] *= CM_IN_KPC;
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
        ("sans-serif", 12).into_font(),
        vec![false],
        None,
        None,
    )
    .unwrap();
}

pub fn core_collapse_background_at_points<T: Fn(f64, f64) -> f64, V: Fn(f64, f64) -> f64>(
    temperature: T,
    neutral_fraction: V,
    rho_s_0: f64,
    r_s_0: f64,
    collapse_progress: f64,
    rho_c: Option<f64>,
    bounds: (f64, f64),
    ang_points: Vec<f64>,
) -> Vec<f64> {
    let rho = parametic_core_collapse(r_s_0, rho_s_0, collapse_progress);

    let r_points = get_r_points(bounds);

    let dark_matter_rho_points = get_rho_points(rho, &r_points);

    let neutral_rho_points = match rho_c {
        None => get_hydrostatic_profile_inwards(&r_points, dark_matter_rho_points, temperature),
        Some(rho_c) => {
            let external_field = get_force_points(dark_matter_rho_points, &r_points);
            get_hydrostatic_profile_outwards(
                &r_points,
                external_field,
                &temperature,
                &neutral_fraction,
                rho_c,
            )
            .0
        }
    };

    let number_density = {
        let mut num_density = Vec::with_capacity(neutral_rho_points.len());
        for i in 0..neutral_rho_points.len() {
            num_density.push(neutral_rho_points[i] / M_PROTON)
        }
        num_density
    };

    let points = {
        let mut points = ang_points.clone();
        for i in 0..points.len() {
            points[i] *= ARC_MIN; // radians
            points[i] *= DISTANCE; // kpc
        }
        points
    };

    let column_density = {
        let mut col_dens = get_column_density_at_points(number_density, &r_points, &points);
        for i in 0..col_dens.len() {
            col_dens[i] /= CM_IN_KPC.powi(2);
        }
        col_dens
    };

    column_density
}

pub fn core_collapse_background<T: Fn(f64, f64) -> f64, V: Fn(f64, f64) -> f64>(
    temperature: T,
    neutral_fraction: V,
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

    let neutral_rho_points = match rho_c {
        None => get_hydrostatic_profile_inwards(&r_points, dark_matter_rho_points, temperature),
        Some(rho_c) => {
            let external_field = get_force_points(dark_matter_rho_points, &r_points);
            get_hydrostatic_profile_outwards(
                &r_points,
                external_field,
                &temperature,
                &neutral_fraction,
                rho_c,
            )
            .0
        }
    };

    //dbg!(&rho_points);

    let number_density = {
        let mut num_density = Vec::with_capacity(neutral_rho_points.len());
        for i in 0..neutral_rho_points.len() {
            let num_dens = neutral_rho_points[i] / M_PROTON; // nH / kpc^-3
            num_density.push(num_dens / CM_IN_KPC.powi(3))
        }
        num_density
    };

    //dbg!(&number_density);

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
            "n_H (num / cm^3)",
            vec![None],
            ("sans-serif", 12).into_font(),
            vec![false],
            None,
            None,
        )
        .unwrap();
    }

    let column_density = {
        let mut col_dens = get_column_density(number_density, &r_points); // nH cm^-3 kpc
        for i in 0..col_dens.len() {
            col_dens[i] *= CM_IN_KPC;
        }
        col_dens
    };

    //dbg!(&column_density);

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
            ("sans-serif", 12).into_font(),
            vec![false],
            None,
            None,
        )
        .unwrap();
    }

    (angular_points, column_density)
}

pub fn parametic_core_collapse(r_s_0: f64, rho_s_0: f64, tau: f64) -> impl Fn(f64) -> f64 {
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

static RELHIC_TN: OnceLock<TnRelation> = OnceLock::new();

pub fn relhic_temperature(_r: f64, rho: f64) -> f64 {
    // https://academic.oup.com/mnras/article/465/4/3913/2544386?login=false Fig 5
    let mut number_density = rho / M_PROTON; // nH kpc^-3
    number_density /= CM_IN_KPC.powi(3); // nH cm^-3

    RELHIC_TN
        .get_or_init(|| TnRelation::from_csv("relhic_Tn_relation.csv").unwrap())
        .temperature_k(number_density)
}

pub fn relhic_neutral_fraction(rho: f64, temperature: f64) -> f64 {
    if temperature < 1e-8 || rho < 1e-8 {
        // this should only happen when rho ~ 0, eta ~0. Really doesn't matter just avoiding 1/T and 1/nH NaNs
        return 0.0;
    }

    // nh must be in cm^-3
    let nh = rho / (M_PROTON * CM_IN_KPC.powi(3));
    // https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2427R/abstract Appendix A2
    let lambda = 315614.0 / temperature;
    let alpha_a = 1.269e-13 * lambda.powf(1.503) / (1.0 + (lambda / 0.522).powf(0.47)).powf(1.923);
    let cap_lambda_t = 1.17e-10 * temperature.sqrt() * (-157809.0 / temperature).exp()
        / (1.0 + (temperature / 1e5).sqrt());
    let gamma_uvb = 2.27e-14; // table 2 line 2
    let gamma_phot = photo_ionization_rate(nh, gamma_uvb);

    let a = alpha_a + cap_lambda_t;
    let b = a + alpha_a + (gamma_phot / nh);
    let c = alpha_a;

    let eta = (b - (b * b - 4.0 * a * c).max(0.0).sqrt()) / (2.0 * a);
    if eta > 1.0 || eta < 0.0 {
        panic!("Nonsense! eta = {eta}")
    }

    eta
}

fn photo_ionization_rate(nh: f64, gamma_uvb: f64) -> f64 {
    // https://ui.adsabs.harvard.edu/abs/2013MNRAS.430.2427R/abstract Appendix A1
    let n_0 = 10_f64.powf(-2.56);
    let alpha_1 = -1.86;
    let alpha_2 = -0.51;
    let beta = 2.83;
    let f = 1.0 - 0.99;

    let ratio = (1.0 - f) * (1.0 + (nh / n_0).powf(beta)).powf(alpha_1)
        + f * (1.0 + (nh / n_0)).powf(alpha_2);

    ratio * gamma_uvb
}

fn get_profile_and_mass<T: Fn(f64, f64) -> f64, V: Fn(f64, f64) -> f64>(
    params: &Vec<f64>,
    r_points: &Vec<f64>,
    temperature: &T,
    neutral_fraction: &V,
) -> (Vec<f64>, f64, Vec<f64>) {
    let dm_rho = parametic_core_collapse(params[1], params[0], params[2]);
    let dm_rho_points = get_rho_points(dm_rho, r_points);
    let external_field = get_force_points(dm_rho_points, r_points);
    let (neutral_rho_gas_points, full_rho_gas_points) = get_hydrostatic_profile_outwards(
        r_points,
        external_field,
        temperature,
        neutral_fraction,
        params[3],
    );
    let total_mass = get_spacial_integral(&full_rho_gas_points, r_points);

    let number_density = {
        let mut num_density = Vec::with_capacity(neutral_rho_gas_points.len());
        for i in 0..neutral_rho_gas_points.len() {
            num_density.push(neutral_rho_gas_points[i] / M_PROTON)
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

    (column_density, total_mass, full_rho_gas_points)
}

pub fn evolution_profile<T: Fn(f64, f64) -> f64, V: Fn(f64, f64) -> f64>(
    anchor_params: &Vec<f64>,
    temperature: &T,
    neutral_fraction: &V,
    data: &fitting::Data,
    font: FontDesc<'static>,
) {
    let halo = Halo::NFW(anchor_params[0], anchor_params[1]);
    let bounds = (1e-5 * INNER_BOUND, 1e2 * halo.r_crit());
    let r_points = get_r_points(bounds);

    // get anchor fit
    let (anchor_gas_rho_points, total_mass, _) =
        get_profile_and_mass(anchor_params, &r_points, temperature, neutral_fraction);
    assert!(is_stable(&anchor_params, temperature, neutral_fraction));

    // get evolution snapshots

    let tau_snapshots = [0.2, 0.4];
    let mut gas_rho_points_vec = Vec::with_capacity(tau_snapshots.len() + 1);
    let mut legends = Vec::with_capacity(tau_snapshots.len() + 1);
    gas_rho_points_vec.push(anchor_gas_rho_points);
    legends.push(Some(format!(
        "Modern Cloud-9, tau = {:.2}",
        anchor_params[2]
    )));
    let mut snapshot_params = anchor_params.clone();
    let mut snapshot_rho_cs: Vec<f64> = Vec::new();
    for tau in tau_snapshots {
        dbg!(tau);
        snapshot_params[2] = tau;
        // Note I don't reset snapshot rho_c between snapshots as snapshots are likely to be closer to each other than the anchor point
        let (mut snapshot_gas_rho_points, mut snapshot_total_mass, _) =
            get_profile_and_mass(&snapshot_params, &r_points, temperature, neutral_fraction);

        let mut percent_diff = (snapshot_total_mass - total_mass) / total_mass;

        let mut iter = 0.0;
        while percent_diff.abs() >= 1e-3 {
            iter += 0.1;
            if iter > 2000.0 {
                panic!("failing to converge for tau={tau}")
            }
            dbg!(percent_diff);
            dbg!(snapshot_params[3]);
            snapshot_params[3] =
                snapshot_params[3].powf((1.0 - percent_diff).powf(1.0 / iter as f64));
            // snapshot_params[3] -= percent_diff * anchor_params[3];

            (snapshot_gas_rho_points, snapshot_total_mass, _) =
                get_profile_and_mass(&snapshot_params, &r_points, temperature, neutral_fraction);

            percent_diff = (snapshot_total_mass - total_mass) / total_mass
        }

        gas_rho_points_vec.push(snapshot_gas_rho_points);
        legends.push(Some(format!("τ = {}", tau)));
        assert!(is_stable(&snapshot_params, temperature, neutral_fraction));
        snapshot_rho_cs.push(snapshot_params[3]);
    }
    dbg!(snapshot_rho_cs);

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
        "Evolution of Cloud-9-like Halo",
        "r (arcmin)",
        "H₁ Column Density (cm⁻²)",
        legends,
        font,
        dashed,
        Some(&data.points),
        Some(&data.y_err),
    )
    .unwrap();
}

pub fn instability_profile<
    T: Fn(f64, f64) -> f64 + std::marker::Sync,
    V: Fn(f64, f64) -> f64 + std::marker::Sync,
>(
    anchor_params: &Vec<f64>,
    bounds: (f64, f64),
    temperature: &T,
    neutral_fraction: &V,
    font: FontDesc<'static>,
) {
    let rho_axis = {
        const AXIS_GRID_NUM: usize = 100;
        let mut rho_axis = Vec::with_capacity(AXIS_GRID_NUM);
        for i in 0..AXIS_GRID_NUM {
            rho_axis.push(
                (bounds.0.ln()
                    + (i as f64) * (bounds.1.ln() - bounds.0.ln()) / ((AXIS_GRID_NUM - 1) as f64))
                    .exp(),
            )
        }
        rho_axis
    };

    let halo = Halo::NFW(anchor_params[0], anchor_params[1]);
    let r_bounds = (1e-5 * INNER_BOUND, 1e3 * halo.r_crit());
    let r_points = get_r_points(r_bounds);

    // get evolution snapshots
    let tau_snapshots = [0.0, 0.3, 0.7, 1.0];
    let mut legends = Vec::with_capacity(tau_snapshots.len() + 1);
    for tau in tau_snapshots {
        legends.push(Some(format!("τ = {}", tau)));
    }
    let mut ms_of_rhos: Vec<Vec<f64>> = tau_snapshots
        .iter()
        .map(|tau| {
            dbg!(tau);
            rho_axis
                .par_iter()
                .map(|rho_c: &f64| {
                    let snapshot_params = vec![anchor_params[0], anchor_params[1], *tau, *rho_c];

                    get_profile_and_mass(&snapshot_params, &r_points, temperature, neutral_fraction)
                        .1
                })
                .collect()
        })
        .collect();

    // get anchor fit
    let (_, total_mass, _) =
        get_profile_and_mass(anchor_params, &r_points, temperature, neutral_fraction);
    ms_of_rhos.push(vec![total_mass; rho_axis.len()]);
    legends.push(Some(format!("Modern Cloud-9 Gas Mass:{:.2e}", total_mass)));

    let mut dashed = vec![true; ms_of_rhos.len()];
    dashed[ms_of_rhos.len() - 1] = false;

    plot_functions(
        &rho_axis,
        &ms_of_rhos,
        "figures/instability.svg",
        "Insability of Cloud-9-like Halo",
        "Gas Central Density",
        "Total Gas Mas",
        legends,
        font,
        dashed,
        None,
        None,
    )
    .unwrap();

    let sparse_rho_axis = {
        const AXIS_GRID_NUM: usize = 10;
        let mut rho_axis = Vec::with_capacity(AXIS_GRID_NUM);
        for i in 0..AXIS_GRID_NUM {
            rho_axis.push(
                (bounds.0.ln()
                    + (i as f64) * (bounds.1.ln() - bounds.0.ln()) / ((AXIS_GRID_NUM - 1) as f64))
                    .exp(),
            )
        }
        rho_axis
    };

    for tau in tau_snapshots {
        let rho_funcs: Vec<Vec<f64>> = sparse_rho_axis
            .iter()
            .map(|rho_c| {
                let mut params = anchor_params.clone();
                params[2] = tau;
                params[3] = *rho_c;
                get_profile_and_mass(&params, &r_points, temperature, neutral_fraction).2
            })
            .collect();

        plot_functions(
            &r_points,
            &rho_funcs,
            &format!("figures/profiles_tau={tau}.svg"),
            "Total Gas Density",
            "r (kpc)",
            "rho (M_s kpc^-3)",
            sparse_rho_axis
                .iter()
                .map(|rho_c| Some(format!("ρ꜀={rho_c:.2e}")))
                .collect(),
            ("sans-serif", 35).into_font(),
            vec![true; rho_funcs.len()],
            None,
            None,
        )
        .unwrap();
    }
}

pub fn instability_showcase<
    T: Fn(f64, f64) -> f64 + std::marker::Sync,
    V: Fn(f64, f64) -> f64 + std::marker::Sync,
>(
    params: &Vec<Vec<f64>>,
    bounds: (f64, f64),
    temperature: &T,
    neutral_fraction: &V,
    font: FontDesc<'static>,
) {
    let rho_axis = {
        const AXIS_GRID_NUM: usize = 100;
        let mut rho_axis = Vec::with_capacity(AXIS_GRID_NUM);
        for i in 0..AXIS_GRID_NUM {
            rho_axis.push(
                (bounds.0.ln()
                    + (i as f64) * (bounds.1.ln() - bounds.0.ln()) / ((AXIS_GRID_NUM - 1) as f64))
                    .exp(),
            )
        }
        rho_axis
    };

    let r_min = 1e-5 * INNER_BOUND;
    let mut r_max = INNER_BOUND;
    for i in 0..params.len() {
        let halo = Halo::NFW(params[i][0], params[i][1]);
        r_max = r_max.max(1e3 * halo.r_crit());
    }
    let r_points = get_r_points((r_min, r_max));

    let mut ms_of_rhos: Vec<Vec<f64>> = params
        .iter()
        .map(|params| {
            rho_axis
                .par_iter()
                .map(|rho_c: &f64| {
                    let snapshot_params = vec![params[0], params[1], params[2], *rho_c];

                    get_profile_and_mass(&snapshot_params, &r_points, temperature, neutral_fraction)
                        .1
                })
                .collect()
        })
        .collect();

    let points: Vec<(f64, f64)> = params
        .iter()
        .map(|params| {
            (
                params[3],
                get_profile_and_mass(&params, &r_points, temperature, neutral_fraction).1,
            )
        })
        .collect();

    let legends: Vec<Option<String>> = params
        .iter()
        .map(|params| {
            let (m200, c200) = rs_rhos_to_m200_c200(params[1], params[0]);
            Some(format!(
                "M₂₀₀ = {:.2e} C₂₀₀ = {:.2e} τ = {:.2}",
                m200, c200, params[2]
            ))
        })
        .rev()
        .collect();

    let dashed = vec![true; ms_of_rhos.len()];

    plot_functions(
        &rho_axis,
        &ms_of_rhos,
        "figures/instability_comp.svg",
        "Insability in Some Cloud-9-like Halos",
        "Gas Central Density",
        "Total Gas Mas",
        legends,
        font,
        dashed,
        Some(&points),
        None,
    )
    .unwrap();
}

pub fn is_stable<T: Fn(f64, f64) -> f64, V: Fn(f64, f64) -> f64>(
    params: &[f64],
    temperature: &T,
    neutral_fraction: &V,
) -> bool {
    let halo = Halo::NFW(params[0], params[1]);
    let r_bounds = (1e-5 * INNER_BOUND, 1e3 * halo.r_crit());
    let r_points = get_r_points(r_bounds);

    let delta = params[3] * 0.02;

    let dm_rho = parametic_core_collapse(params[1], params[0], params[2]);
    let dm_rho_points = get_rho_points(dm_rho, &r_points);
    let external_field = get_force_points(dm_rho_points, &r_points);

    let (_neutral_rho_gas_points, lower_full_rho_gas_points) = get_hydrostatic_profile_outwards(
        &r_points,
        external_field.clone(),
        temperature,
        neutral_fraction,
        params[3] - delta,
    );
    let lower_total_mass = get_spacial_integral(&lower_full_rho_gas_points, &r_points);
    let (_neutral_rho_gas_points, upper_full_rho_gas_points) = get_hydrostatic_profile_outwards(
        &r_points,
        external_field.clone(),
        temperature,
        neutral_fraction,
        params[3] + delta,
    );
    let upper_total_mass = get_spacial_integral(&upper_full_rho_gas_points, &r_points);

    upper_total_mass > lower_total_mass
}

#[cfg(test)]
mod tests {
    use super::*;
    use plotters::style::FontDesc;
    use plotters::style::{FontFamily, FontStyle};

    #[test]
    fn hydro_test() {
        const TEMP: f64 = 1e4;
        let temperature = |_: f64, _: f64| -> f64 { TEMP };
        let neutral_fraction = |_: f64, _: f64| -> f64 { 0.0 };
        let mu = 0.5;
        let r_bounds = (1e-1 * INNER_BOUND, 1e10 * INNER_BOUND);
        let r_points = get_r_points(r_bounds);

        let bounds: (f64, f64) = (1e0, 1e7);
        let rho_axis = {
            const AXIS_GRID_NUM: usize = 100;
            let mut rho_axis = Vec::with_capacity(AXIS_GRID_NUM);
            for i in 0..AXIS_GRID_NUM {
                rho_axis.push(
                    (bounds.0.ln()
                        + (i as f64) * (bounds.1.ln() - bounds.0.ln())
                            / ((AXIS_GRID_NUM - 1) as f64))
                        .exp(),
                )
            }
            rho_axis
        };

        let rho_funcs: Vec<Vec<f64>> = rho_axis
            .par_iter()
            .map(|rho_c| {
                let (_neutral_rho_gas_points, full_rho_gas_points) =
                    get_hydrostatic_profile_outwards(
                        &r_points,
                        vec![0.0; r_points.len()],
                        temperature,
                        neutral_fraction,
                        *rho_c,
                    );

                full_rho_gas_points
            })
            .collect();

        let masses: Vec<f64> = rho_funcs
            .par_iter()
            .map(|rho_points| get_spacial_integral(rho_points, &r_points))
            .collect();

        let analytical_points: Vec<f64> = r_points
            .iter()
            .map(|r| {
                let rho = (TEMP / (mu * MP_OVER_KB)) / (2.0 * PI * GG * r.powi(2));
                // [kpc^2 s^-2][km^-2 kpc^-1 Msun s^2][kpc^-2] = M_sun km^-2 kpc^-1
                rho * KM_IN_KPC.powi(2)
            })
            .collect();

        let mut plot_funcs = vec![analytical_points.clone()];
        plot_funcs.extend(rho_funcs.clone());
        let mut legends: Vec<Option<String>> = rho_axis
            .iter()
            .map(|rho_c| Some(format!("{rho_c:.2e}")))
            .collect();
        legends.insert(0, Some(format!("lim")));
        let mut dashed = vec![false];
        dashed.extend(vec![true; rho_axis.len()]);
        plot_functions(
            &r_points,
            &plot_funcs,
            "diagnostics/hydro_test.svg",
            "Hydro tests",
            "r (kpc)",
            "rho (M_s kpc^-3)",
            legends,
            ("sans-serif", 35).into_font(),
            dashed,
            None,
            None,
        )
        .unwrap();

        let mut analytical_mass_asymptote =
            2.0 * (TEMP / (mu * MP_OVER_KB)) * r_points[r_points.len() - 1] / GG;
        // [kpc^2 s^-2][km^-2 kpc^-1 Msun s^2][kpc] = kpc^2 Msun km^-2
        analytical_mass_asymptote *= KM_IN_KPC.powi(2);

        plot_functions(
            &rho_axis,
            &vec![
                masses.clone(),
                vec![analytical_mass_asymptote; rho_axis.len()],
            ],
            "diagnostics/M_vs_rho_test.svg",
            "Hydro tests",
            "rho_c (M_s kpc^-3)",
            "M (M_s)",
            vec![Some(format!("Numerical")), Some(format!("Analytic"))],
            ("sans-serif", 35).into_font(),
            vec![false, true],
            None,
            None,
        )
        .unwrap();

        const REACHBACK: usize = 10;
        let asymptote_checks: Vec<f64> = rho_funcs
            .par_iter()
            .map(|rho_points| {
                let delta_log_rho = rho_points[rho_points.len() - 1].ln()
                    - rho_points[rho_points.len() - REACHBACK].ln();
                let delta_log_r =
                    r_points[r_points.len() - 1].ln() - r_points[r_points.len() - REACHBACK].ln();
                delta_log_rho / delta_log_r
            })
            .collect();

        for log_slope in asymptote_checks {
            if (log_slope + 2.0).abs() > 0.01 {
                dbg!(log_slope);
                panic!("Log slope not correctly on asymptote!")
            }
        }

        let square_err: f64 = rho_funcs
            .iter()
            .map(|rho_points| {
                (rho_points[rho_points.len() - REACHBACK]
                    - analytical_points[analytical_points.len() - REACHBACK])
                    / analytical_points[analytical_points.len() - REACHBACK]
            })
            .map(|err| err.powi(2))
            .sum();

        let rms = (square_err / rho_axis.len() as f64).sqrt();
        if rms > 0.01 {
            panic!("Large radii asymptotes not quite right: rms from analytical: {rms}")
        }

        let square_mass_err: f64 = masses
            .iter()
            .map(|mass| (mass - analytical_mass_asymptote) / analytical_mass_asymptote)
            .map(|err| err.powi(2))
            .sum();

        let mass_rms = (square_mass_err / rho_axis.len() as f64).sqrt();
        if mass_rms > 0.01 {
            panic!("Mass asymptote not matching analytical solution")
        }
    }

    #[test]
    fn harmonic_test() {
        const TEMP: f64 = 1e10;
        let temperature = |_: f64, _: f64| -> f64 { TEMP };
        let neutral_fraction = |_: f64, _: f64| -> f64 { 0.0 };
        let r_bounds = (1e-3 * INNER_BOUND, 1e10 * INNER_BOUND);
        let r_points = get_r_points(r_bounds);

        const OMEGA: f64 = 1e-15;

        let external_field: Vec<f64> = r_points.iter().map(|r| -OMEGA.powi(2) * r).collect();

        // start at low rho_c and high T to effectively turn off self gravity
        const RHO_C: f64 = 1.0;
        let (_, rho_points) = get_hydrostatic_profile_outwards(
            &r_points,
            external_field,
            temperature,
            neutral_fraction,
            RHO_C,
        );

        let analytic_points: Vec<f64> = r_points
            .iter()
            .map(|r| RHO_C * (-(0.5 * MP_OVER_KB * OMEGA.powi(2) * r.powi(2)) / (2.0 * TEMP)).exp())
            .collect();

        plot_functions(
            &r_points,
            &vec![rho_points.clone(), analytic_points.clone()],
            "diagnostics/harmonic_test.svg",
            "Hydro tests",
            "r (kpc)",
            "rho (M_s kpc^-3)",
            vec![Some(format!("Numerical")), Some(format!("Analytic"))],
            ("sans-serif", 35).into_font(),
            vec![false, true],
            None,
            None,
        )
        .unwrap();

        let square_err: f64 = rho_points
            .iter()
            .zip(analytic_points)
            .map(|(numerical, analytic)| (numerical - analytic).powi(2))
            .sum();

        let rms = (square_err / rho_points.len() as f64).sqrt();

        if rms > 0.001 {
            panic!("RMS error too high, something wrong in hydrostatic integrator. RMS: {rms}")
        }
    }

    #[test]
    fn test_spacial_integral() {
        let r_bounds = (1e-3 * INNER_BOUND, 1e10 * INNER_BOUND);
        let r_points = get_r_points(r_bounds);

        let numerical = get_spacial_integral(&vec![1.0; r_points.len()], &r_points);
        let analytical = 4.0 * PI * r_bounds.1.powi(3) / 3.0;

        let err = (numerical - analytical) / analytical;

        if err.abs() > 1e-12 {
            panic!("Spacial integral isn't matching up: %err = {err}");
        }
    }

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
                ("sans-serif", 12).into_font(),
                vec![false],
                None,
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
                ("sans-serif", 12).into_font(),
                vec![false],
                None,
                None,
            )
            .unwrap();
            panic!("rms error too high! rms_err = {rms_err}")
        }
    }

    #[test]
    fn relhic_neutral_frac_plot() {
        let font = FontDesc::new(FontFamily::SansSerif, 16.0, FontStyle::Normal);

        let mut t_points = Vec::new();
        let mut neut_frac_points_t = Vec::new();
        let ref_n = 2e-3;
        let ref_rho = (M_PROTON * CM_IN_KPC.powi(3)) * ref_n;

        let (t_min, t_max): (f64, f64) = (1e3, 1e5);
        for i in 0..1000 {
            let progress = i as f64 / 999.0;
            t_points.push((t_min.ln() + (progress) * (t_max.ln() - t_min.ln())).exp());
            neut_frac_points_t.push(relhic_neutral_fraction(ref_rho, t_points[i]));
        }

        plot_functions(
            &t_points,
            &vec![neut_frac_points_t],
            "T_vs_neutral_frac.svg",
            &format!("T_vs_neutral_frac at rho={ref_rho}"),
            "T",
            "Neutral Frac",
            vec![None],
            font.clone(),
            vec![true],
            None,
            None,
        )
        .unwrap();

        let mut rho_points = Vec::new();
        let mut neut_frac_points_rho = Vec::new();
        let ref_temp = 1e4;

        let (rho_min, rho_max): (f64, f64) = (5e3, 2e5);
        for i in 0..1000 {
            let progress = i as f64 / 999.0;
            rho_points.push((rho_min.ln() + (progress) * (rho_max.ln() - rho_min.ln())).exp());
            neut_frac_points_rho.push(relhic_neutral_fraction(rho_points[i], ref_temp));
        }

        plot_functions(
            &rho_points,
            &vec![neut_frac_points_rho],
            "rho_vs_neutral_frac.svg",
            &format!("rho_vs_neutral_frac at T={ref_temp}"),
            "rho",
            "Neutral Frac",
            vec![None],
            font,
            vec![true],
            None,
            None,
        )
        .unwrap();
    }

    #[test]
    fn core_collapse_reproduce_nfw() {
        let (rho_s_0, r_s_0) = (1e6, 4.2);
        let tau: f64 = 0.0;
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

        dbg!(rho_s, r_s, r_c);

        if rho_s - rho_s_0 > 1e-5 || r_s - r_s_0 > 1e-5 || r_c > 1e-5 {
            panic!("cored profile does not recover nfw at tau = 0")
        }
    }
}
