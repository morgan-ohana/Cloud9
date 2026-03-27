use plotters::style::IntoFont;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;
use rayon::prelude::*;

use crate::halo::{Halo, McrSource, m200_c200_to_rs_rhos, mass_concentration_relation};
use crate::hydrostatics::{isothermal_core_collapse_background_at_points, relhic_temperature};
use crate::{hydrostatics::isothermal_core_collapse_background, plotting::plot_functions};

#[derive(Clone)]
struct Walker {
    params: [f64; 4],
    log_prob: f64,
    rng: Pcg64,
}

impl Walker {
    fn lin_params(&self) -> [f64; 4] {
        let mut lin_params = self.params.clone();
        for s in 0..4 {
            if LOG_STEP[s] {
                lin_params[s] = 10.0_f64.powf(self.params[s])
            }
        }

        lin_params
    }
}

fn sample_z(rng: &mut impl Rng, a: f64) -> f64 {
    // pdf ~ 1/sqrt(z) [1/a, a]
    // cdf ~ 2sqrt(z) + C [1/a, a]
    // Normalizing on the range it must be: cdf = (sqrt(z) - 1/sqrt(a)) / (sqrt(a) - 1/sqrt(a))
    // u is random sample from cdf, so u * (sqrt(a) - 1/sqrt(a)) + 1/sqrt(a) = sqrt(z)
    // => z = (1/sqrt(a) + u * (sqrt(a) - 1/sqrt(a)))^2
    let u: f64 = rng.random();
    let b = 1.0 / a.sqrt();
    (b + u * (a.sqrt() - b)).powi(2)
}

fn stretch_move_parallel(
    walkers: &mut [Walker],
    log_likelihood: &(impl Fn([f64; 4]) -> f64 + std::marker::Sync),
    a: f64,
    bounds: &[[f64; 2]; 4],
) {
    let n = walkers.len();
    let d = 4;

    // Copy walkers so reads are consistent during parallel update
    let walkers_old = walkers.to_vec();

    walkers.par_iter_mut().enumerate().for_each(|(i, walker)| {
        // choose partner
        let mut j = walker.rng.random_range(0..n);
        while j == i {
            j = walker.rng.random_range(0..n);
        }

        let z = sample_z(&mut walker.rng, a);
        let mut proposal = walker.params;

        for k in 0..d {
            proposal[k] =
                walkers_old[j].params[k] + z * (walker.params[k] - walkers_old[j].params[k]);

            if proposal[k] < bounds[k][0] || proposal[k] > bounds[k][1] {
                // if out of bounds just stop.
                // dbg!("FUCK");
                // dbg!(k);
                // Walker won't move so this will correctly count as a rejection
                return;
            }
        }

        let proposed_log_prob = log_likelihood(proposal);

        let log_accept = (d as f64 - 1.0) * z.ln() + proposed_log_prob - walker.log_prob;

        if log_accept >= 0.0 || walker.rng.random::<f64>() < log_accept.exp() {
            walker.params = proposal;
            walker.log_prob = proposed_log_prob;
        }
    });
}

const LOG_STEP: [bool; 4] = [false, false, false, false];
const LOG_PRIOR: [bool; 4] = [false, false, false, false];
pub fn find_parameters_mcmc(
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
    initial_guess: [f64; 4],
    lin_bounds: &[[f64; 2]; 4],
    num_steps: usize,
    burn_in: usize,
    n_walkers: usize,
    prior: Prior,
) -> ([f64; 4], Vec<[f64; 4]>, Vec<f64>) {
    let bounds = {
        let mut bounds = lin_bounds.clone();
        for s in 0..4 {
            if LOG_STEP[s] {
                bounds[s][0] = bounds[s][0].log10();
                bounds[s][1] = bounds[s][1].log10();
            }
        }

        bounds
    };

    let mut rng = Pcg64::seed_from_u64(42);
    let a = 2.0;

    // initialize walkers
    let mut walkers = Vec::with_capacity(n_walkers);
    for i in 0..n_walkers {
        let mut p = initial_guess;
        for s in 0..4 {
            if LOG_STEP[s] {
                p[s] = p[s].log10()
            }
            p[s] *= 1.0 + 0.1 * (rng.random::<f64>() - 0.5);
        }
        let lp = log_likelihood(p, data, y_error_bar, true, prior.clone());
        walkers.push(Walker {
            params: p,
            log_prob: lp,
            rng: Pcg64::seed_from_u64(42 + i as u64 * 7919),
        });
    }

    let mut chains = vec![Vec::new(); n_walkers];
    let mut likelihoods = vec![Vec::new(); n_walkers];

    for step in 0..num_steps {
        stretch_move_parallel(
            &mut walkers,
            &|p| log_likelihood(p, data, y_error_bar, true, prior.clone()),
            a,
            &bounds,
        );

        if step >= burn_in {
            for w in 0..n_walkers {
                chains[w].push(walkers[w].lin_params());
                likelihoods[w].push(walkers[w].log_prob);
            }
        }

        if step % 1000 == 0 {
            println!("Step {}", step);
        }
    }

    let (combined_chain, combined_likelihoods) = combine_chains((chains, likelihoods));

    // Best params
    let best_idx = combined_likelihoods
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let best_params = combined_chain[best_idx];

    (best_params, combined_chain, combined_likelihoods)
}

pub fn combine_chains(
    (chains, likelihoods): (Vec<Vec<[f64; 4]>>, Vec<Vec<f64>>),
) -> (Vec<[f64; 4]>, Vec<f64>) {
    let mut combined_chain = Vec::new();
    let mut combined_likelihoods = Vec::new();

    for c in 0..chains.len() {
        combined_chain.extend_from_slice(&chains[c]);
        combined_likelihoods.extend_from_slice(&likelihoods[c]);
    }

    (combined_chain, combined_likelihoods)
}

pub fn split_chains(
    chain: &([f64; 4], Vec<[f64; 4]>, Vec<f64>),
    length: usize,
) -> Vec<([f64; 4], Vec<[f64; 4]>, Vec<f64>)> {
    let num_chains = chain.1.len() / length;

    // Ensure the combined data length is consistent
    assert_eq!(chain.1.len() % length, 0);
    assert_eq!(chain.2.len(), chain.1.len());

    // Prepare containers for the split chains
    let mut split_params = vec![Vec::with_capacity(length); num_chains];
    let mut split_likelihoods = vec![Vec::with_capacity(length); num_chains];

    // Distribute the combined parameter vectors and likelihoods
    for (i, params) in chain.1.chunks(length).enumerate() {
        split_params[i].extend_from_slice(params);
    }
    for (i, liks) in chain.2.chunks(length).enumerate() {
        split_likelihoods[i].extend_from_slice(liks);
    }

    // Build the output chains, each with its own best parameters
    let mut result = Vec::with_capacity(num_chains);
    for i in 0..num_chains {
        // Find index of maximum likelihood in this chain
        let best_idx = split_likelihoods[i]
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(idx, _)| idx)
            .expect("Chain should not be empty");

        let best_params = split_params[i][best_idx];
        result.push((
            best_params,
            std::mem::take(&mut split_params[i]),
            std::mem::take(&mut split_likelihoods[i]),
        ));
    }

    result
}

pub fn calculate_gelman_rubin(chains: &[([f64; 4], Vec<[f64; 4]>, Vec<f64>)]) -> Vec<f64> {
    let m = chains.len() as f64;
    let n = chains[0].1.len() as f64; // Samples per chain after burn-in

    let mut r_hat = vec![0.0; 4];

    for param_idx in 0..4 {
        // Collect samples for this parameter
        let mut param_samples = Vec::new();
        for (_, chain, _) in chains {
            let chain_samples: Vec<f64> = chain.iter().map(|p| p[param_idx]).collect();
            param_samples.push(chain_samples);
        }

        // Calculate within-chain variance
        let mut chain_means = Vec::new();
        let mut within_var = 0.0;

        for samples in &param_samples {
            let mean: f64 = samples.iter().sum::<f64>() / n;
            chain_means.push(mean);

            let variance: f64 =
                samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
            within_var += variance;
        }
        within_var /= m;

        // Calculate between-chain variance
        let overall_mean: f64 = chain_means.iter().sum::<f64>() / m;

        let between_var: f64 = chain_means
            .iter()
            .map(|&mean| (mean - overall_mean).powi(2))
            .sum::<f64>()
            / (m - 1.0);

        // Calculate pooled variance
        let pooled_var = ((n - 1.0) / n) * within_var + between_var;

        // R-hat statistic (should approach 1.0 as chains converge)
        r_hat[param_idx] = (pooled_var / within_var).sqrt();
    }

    r_hat
}

const NAMES: [&str; 4] = ["m200", "c200", "tau", "rho_c"];
pub fn likelihood_slice_profile(
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
    anchor_point: [f64; 4],
    slice_idx: usize,
    bounds: [f64; 2],
    prior: Prior,
) {
    const N: usize = 100;
    let mut x_range = Vec::with_capacity(N);
    let mut likelihoods = Vec::with_capacity(N);
    for n in 0..N {
        let x = match LOG_STEP[slice_idx] {
            true => {
                (bounds[0].ln() + (n as f64 / N as f64) * (bounds[1].ln() - bounds[0].ln())).exp()
            }
            false => bounds[0] + (n as f64 / N as f64) * (bounds[1] - bounds[0]),
        };

        let mut params = anchor_point.clone();
        params[slice_idx] = x;
        for s in 0..4 {
            if LOG_STEP[s] {
                params[s] = params[s].log10();
            }
        }

        let likelihood = log_likelihood(params, data, y_error_bar, true, prior.clone());

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

#[derive(Clone, Debug)]
pub enum Prior {
    MassConcentrationRelation(McrSource),
    None,
}

fn log_likelihood(
    mut params: [f64; 4],
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
    m200_input: bool,
    prior: Prior,
) -> f64 {
    // undo log scale and compute jacobian
    let mut log_jacobian = 0.0;
    for s in 0..4 {
        if LOG_STEP[s] {
            params[s] = 10.0_f64.powf(params[s]);

            if !LOG_PRIOR[s] {
                log_jacobian += params[s].ln();
            }
        }
    }

    let log_prior_likelihood = match prior {
        Prior::MassConcentrationRelation(mcr_source) => {
            let mut mcr_prior = 0.0;
            if m200_input {
                // received as m200, c200, tau, rho_c

                // compute prior biases from mass-concetration-relation
                let (mean_log10c, sigma_log10c) =
                    mass_concentration_relation(params[0], mcr_source);
                let log10c = params[1].log10();
                mcr_prior = -0.5 * ((log10c - mean_log10c) / sigma_log10c).powi(2);
            } else {
                // No implemented Mcr prior for rho_s r_s space.
            }

            mcr_prior
        }
        Prior::None => 0.0,
    };

    if m200_input {
        // recieved as m200, c200, tau, rho_c
        (params[1], params[0]) = m200_c200_to_rs_rhos(params[0], params[1]);
    }

    let halo = Halo::NFW(params[0], params[1]);

    let ang_points = data.iter().map(|(ang, _)| *ang).collect();

    // must be passed as rho_s_0, r_s_0, tau, rho_c
    let model_points = isothermal_core_collapse_background_at_points(
        relhic_temperature,
        params[0],
        params[1],
        params[2],
        Some(params[3]),
        (0.1 * data[0].0, halo.r_crit()), // Must go out to r_crit so I have whole halo for projected density
        ang_points,
    );

    // compute X^2
    let mut chi_squared = 0.0;
    let mut sum_ln_sigma = 0.0;
    for i in 0..data.len() {
        let diff = model_points[i] - data[i].1;
        let spread = match diff > 0.0 {
            true => y_error_bar[i].1 - data[i].1,
            false => data[i].1 - y_error_bar[i].0,
        };
        chi_squared += (diff / spread).powi(2);
        sum_ln_sigma += spread.ln();
    }

    -0.5 * chi_squared - sum_ln_sigma + log_prior_likelihood + log_jacobian
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
    let fit = isothermal_core_collapse_background(
        relhic_temperature,
        params[0],
        params[1],
        params[2],
        rho_c,
        (0.1 * data[0].0, halo.r_crit()),
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
