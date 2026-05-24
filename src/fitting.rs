use ensemble_mcmc::*;
use plotters::style::IntoFont;
use std::f64::consts::PI;

use crate::constants::*;
use crate::halo::{
    Halo, McrSource, m200_c200_to_rs_rhos, mass_concentration_relation, rs_rhos_to_m200_c200,
};
use crate::hydrostatics::{
    core_collapse_background_at_points, relhic_neutral_fraction, relhic_temperature,
};
// use crate::mcmc::*;
use crate::{hydrostatics::core_collapse_background, plotting::plot_functions};

#[derive(Clone)]
pub struct Data {
    pub points: Vec<(f64, f64)>,
    pub y_err: Vec<(f64, f64)>,
}

impl Data {
    pub fn init() -> Self {
        // Digitized from https://doi.org/10.3847/1538-4357/ad65d9 figure 4

        let points: Vec<(f64, f64)> = vec![
            (0.15649748494408605, 15595734576138530000.0),
            (0.18722388772639217, 15290963190523656000.0),
            (0.2252243432550796, 14924613361774703000.0),
            (0.26936457583283363, 14424784086219026000.0),
            (0.32301280500582225, 13894659794863262000.0),
            (0.38799031821006913, 13091285860285037000.0),
            (0.4647330501352986, 12151824423121148000.0),
            (0.5571917155289376, 11103690653134041000.0),
            (0.6692414165028057, 9918955303017603000.0),
            (0.801419168108796, 8606124456897346000.0),
            (0.9595907722097119, 7014426686842890000.0),
            (1.1500508076747187, 5250966932011177000.0),
            (1.3771319397378698, 3511010519320841000.0),
            (1.6463993031339876, 2293590617222256000.0),
            (1.983419488324015, 1615272370298475800.0),
        ];

        let mut y_err: Vec<(f64, f64)> = vec![
            (13732078146520720000.0, 17466755309837742000.0),
            (13427854263419597000.0, 17229242059185019000.0),
            (13003516161985374000.0, 16866179741318795000.0),
            (12461619020702566000.0, 16448095066290774000.0),
            (11695185329352206000.0, 15955195009096049000.0),
            (10716711258649842000.0, 15359790383829522000.0),
            (9581053603324199000.0, 14708083512691229000.0),
            (8205955659163100000.0, 13903084887830374000.0),
            (6793221053775176000.0, 12983171214351038000.0),
            (5337519460514788000.0, 11805617040307792000.0),
            (3880856839277999600.0, 10151147099502905000.0),
            (2370604088405835000.0, 8104671612570444000.0),
            (684370527200647700.0, 6345004385610885000.0),
            (6e17, 5195521117287574000.0), // lower bound here doesn't matter, will be manually set below since not visible on plot
            (6e17, 4312913465944557000.0), // ditto
        ];

        // Lower bound for unknown ones taken symmetrical in linear space
        let n = y_err.len();
        y_err[n - 1].0 = 2.0 * points[n - 1].1 - y_err[n - 1].1;
        y_err[n - 2].0 = 2.0 * points[n - 2].1 - y_err[n - 2].1;

        Self { points, y_err }
    }
}
const NAMES: [&str; 4] = ["m200", "c200", "tau", "rho_c"];
pub fn likelihood_slice_profile(
    data: &Data,
    anchor_point: [f64; 4],
    slice_idx: usize,
    bounds: [f64; 2],
    prior: Prior,
) {
    const N: usize = 100;
    let mut x_range = Vec::with_capacity(N);
    let mut likelihoods = Vec::with_capacity(N);
    for n in 0..N {
        let x = bounds[0] + (n as f64 / N as f64) * (bounds[1] - bounds[0]);

        let mut params = anchor_point.clone();
        params[slice_idx] = x;

        let likelihood = log_likelihood_full(&params[..], data, true, prior.clone());

        x_range.push(x);
        likelihoods.push(likelihood);
    }

    plot_functions(
        &x_range,
        &vec![likelihoods],
        "likelihood_profile.png",
        "Likelihood Profile",
        NAMES[slice_idx],
        "log likelihood",
        vec![None],
        ("sans-serif", 12).into_font(),
        vec![false],
        None,
        None,
    )
    .unwrap()
}

#[derive(Clone, Debug)]
pub enum Prior {
    MassConcentrationRelation(McrSource),
    None,
}

pub enum Cloud9MCMCCore {
    Full(Data, Prior, [[f64; 2]; 4]),
    FixedCrossSection(Data, Prior, [[f64; 2]; 3], f64),
}

impl Cloud9MCMCCore {
    pub fn init(
        data: Data,
        prior: Prior,
        bounds: [[f64; 2]; 4],
        fixed_cross_section: Option<f64>,
    ) -> Self {
        match fixed_cross_section {
            None => Cloud9MCMCCore::Full(data, prior, bounds),
            Some(sigma) => {
                let no_tau_bounds = [bounds[0], bounds[1], bounds[3]];
                Cloud9MCMCCore::FixedCrossSection(data, prior, no_tau_bounds, sigma)
            }
        }
    }
}

impl MCMCCore for Cloud9MCMCCore {
    fn get_bounds(&self) -> &[[f64; 2]] {
        match self {
            Cloud9MCMCCore::Full(_, _, bounds) => &bounds[..],
            Cloud9MCMCCore::FixedCrossSection(_, _, bounds, _sigma) => &bounds[..],
        }
    }

    fn get_log_likelihood(&self, params: &[f64]) -> f64 {
        match self {
            Cloud9MCMCCore::Full(data, prior, _) => {
                log_likelihood_full(params, data, true, prior.clone())
            }
            Cloud9MCMCCore::FixedCrossSection(data, prior, _, cross_section) => {
                assert_eq!(params.len(), 3);

                const AGE: f64 = 10.0;
                let (r_s, rho_s) = m200_c200_to_rs_rhos(params[0], params[1]);

                // cm^2 g^-1 Gyr kpc Ms kpc^-3 (km^2 kpc Ms^-1 s^-2 Ms kpc^-3)^0.5 = cm^2 g^-1 Gyr Ms kpc^-2 km kpc^-1 s^-1 = (cm/kpc)^2 (km/kpc) (Gyr/s) (Ms/g)
                let tau =
                    (0.75 * cross_section * AGE * r_s * rho_s * (4.0 * PI * GG * rho_s).sqrt()
                        / 150.0)
                        * S_IN_GYR
                        * G_IN_MSUN
                        / (KM_IN_KPC * CM_IN_KPC.powi(2));

                // if tau > 1.0 {
                //     println!(
                //         "tau = {tau} with M200 = {}, C200 = {}",
                //         params[0], params[1]
                //     )
                // }

                log_likelihood_full(&[rho_s, r_s, tau, params[2]], data, false, prior.clone())
            }
        }
    }
}

fn log_likelihood_full(params: &[f64], data: &Data, m200_input: bool, prior: Prior) -> f64 {
    assert_eq!(params.len(), 4);

    let log_prior_likelihood = match prior {
        Prior::MassConcentrationRelation(mcr_source) => {
            if m200_input {
                // received as m200, c200, tau, rho_c
                let (mean_log10c, sigma_log10c) =
                    mass_concentration_relation(params[0], mcr_source);
                let log10c = params[1].log10();
                -0.5 * ((log10c - mean_log10c) / sigma_log10c).powi(2)
            } else {
                // recieved as rhos rs tau rho_c
                let (m200, c200) = rs_rhos_to_m200_c200(params[1], params[0]);
                let (mean_log10c, sigma_log10c) = mass_concentration_relation(m200, mcr_source);
                let log10c = c200.log10();
                -0.5 * ((log10c - mean_log10c) / sigma_log10c).powi(2)
            }
        }
        Prior::None => 0.0,
    };

    let (r_s, rho_s) = if m200_input {
        // recieved as m200, c200, tau, rho_c
        m200_c200_to_rs_rhos(params[0], params[1])
    } else {
        (params[1], params[0])
    };

    let halo = Halo::NFW(rho_s, r_s);

    let ang_points = data.points.iter().map(|(ang, _)| *ang).collect();

    // must be passed as rho_s_0, r_s_0, tau, rho_c
    let model_points = core_collapse_background_at_points(
        relhic_temperature,
        relhic_neutral_fraction,
        rho_s,
        r_s,
        params[2],
        Some(params[3]),
        (INNER_BOUND, halo.r_crit()), // Must go out to r_crit so I have whole halo for projected density
        ang_points,
    );

    // compute X^2
    let mut chi_squared = 0.0;
    let mut sum_ln_sigma = 0.0;
    for i in 0..data.points.len() {
        let diff = model_points[i] - data.points[i].1;
        let spread = match diff > 0.0 {
            true => data.y_err[i].1 - data.points[i].1,
            false => data.points[i].1 - data.y_err[i].0,
        };
        chi_squared += (diff / spread).powi(2);
        sum_ln_sigma += spread.ln();
    }

    -0.5 * chi_squared - sum_ln_sigma + log_prior_likelihood
}

pub fn find_parameters_gradient_descent(
    data: &Vec<(f64, f64)>,
    initial_guess: [f64; 4],
    fix_tau: Option<f64>,
    inwards: bool,
) -> [f64; 4] {
    let mut params = initial_guess;

    let mut error = Vec::new();
    error.push(1.0);
    error.push(1.0);
    let mut rel_error_change = f64::MAX;
    let mut n = 0;

    while n < 10 || (error[error.len() - 1] > 0.01 && rel_error_change > 0.0 && n < 2000) {
        error.push(get_rms_err_of_fit(params, &data, false, inwards));
        rel_error_change =
            ((error[error.len() - 1] - error[error.len() - 2]) / error[error.len() - 2]).abs();

        let mut grad: [f64; 4] = [0.0; 4];

        let grad_scale = (0.1 * (1.0_f64 - 5e-3).powi(n)).max(1e-3);
        for i in 0..params.len() {
            let mut increased_params = params.clone();
            increased_params[i] += 0.5 * grad_scale * params[i];
            let increase_error = get_rms_err_of_fit(increased_params, &data, false, inwards);

            let mut decreased_params = params.clone();
            decreased_params[i] -= 0.5 * grad_scale * params[i];
            let decrease_error = get_rms_err_of_fit(decreased_params, &data, false, inwards);
            grad[i] = (increase_error - decrease_error) / (grad_scale * params[i])
        }

        for i in 0..params.len() {
            params[i] -= grad[i] * params[i].powi(2) * grad_scale
        }

        if let Some(fixed_tau) = fix_tau {
            params[2] = fixed_tau;
        }

        dbg!(n);
        //dbg!(grad);
        //dbg!(params);
        dbg!(error[error.len() - 1]);
        dbg!(rel_error_change);
        n += 1
    }

    plot_functions(
        &(0..n).map(|n| n as f64).collect(),
        &vec![error],
        "fit_error.png",
        "Fitting Error",
        "Steps",
        "Error",
        vec![None],
        ("sans-serif", 12).into_font(),
        vec![false],
        None,
        None,
    )
    .unwrap();

    params
}

fn get_rms_err_of_fit(
    mut params: [f64; 4],
    data: &Vec<(f64, f64)>,
    m200_input: bool,
    inwards: bool,
) -> f64 {
    if m200_input {
        // recieved as m200, c200, tau, rho_c
        (params[1], params[0]) = m200_c200_to_rs_rhos(params[0], params[1]);
    }

    let halo = Halo::NFW(params[0], params[1]);

    let rho_c = match inwards {
        true => None,
        false => Some(params[3]),
    };

    // must be passed as rho_s_0, r_s_0, tau, rho_c
    let fit = core_collapse_background(
        relhic_temperature,
        relhic_neutral_fraction,
        params[0],
        params[1],
        params[2],
        rho_c,
        (INNER_BOUND, halo.r_crit()),
        false,
    );

    let mut error = 0.0;
    for point in data {
        let mut high = fit.0.len() - 1;
        let mut low = 0;
        while high - low > 1 {
            let mid = (low + high) / 2;
            if point.0 > fit.0[mid] {
                low = mid
            } else {
                high = mid
            }
        }

        error += ((fit.1[low] - point.1) / point.1).powi(2);
    }

    error /= data.len() as f64;
    error.sqrt()
}

pub fn calculate_statistics(chain: &Vec<[f64; 4]>, best_params: &[f64; 4]) -> [f64; 4] {
    let n = chain.len() as f64;

    println!("\nParameter Statistics:");
    println!(
        "{:<10} {:<12} {:<12} {:<12}",
        "Param", "Best", "Mean", "Std Dev"
    );
    println!("{}", "-".repeat(48));

    let names = ["m200_0", "c200_0", "tau", "rho_c"];
    let mut mean = [0.0; 4];

    for i in 0..4 {
        mean[i] = chain.iter().map(|p| p[i]).sum::<f64>() / n;
        let variance = chain.iter().map(|p| (p[i] - mean[i]).powi(2)).sum::<f64>() / n;
        let std_dev = variance.sqrt();

        println!(
            "{:<10} {:<12} {:<12} {:<12}",
            names[i],
            format_number(best_params[i]),
            format_number(mean[i]),
            format_number(std_dev)
        );
    }

    mean
}

fn format_number(x: f64) -> String {
    // Use scientific notation for very large or very small numbers
    if x.abs() >= 1e3 || (x.abs() > 0.0 && x.abs() < 1e-2) {
        format!("{:.4e}", x)
    } else {
        format!("{:.4}", x)
    }
}
