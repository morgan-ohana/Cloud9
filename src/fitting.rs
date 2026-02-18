use core::num;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};
use rand_pcg::Pcg64;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::constants::UVB_TEMP;
use crate::halo::{Halo, m200_c200_to_rs_rhos};
use crate::plotting::create_chain_trace_plots;
use crate::{hydrostatics::isothermal_core_collapse_background, plotting::plot_function};

pub fn find_parameters_mcmc(
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
    initial_guess: [f64; 4],
    bounds: [(f64, f64); 4],
    fix_tau: Option<f64>,
    num_steps: usize,
    burn_in: usize,
    inwards: bool,
    num_chains: usize,
) -> Result<([f64; 4], Vec<[f64; 4]>, Vec<f64>), Box<dyn std::error::Error>> {
    // Run parallel chains
    let (chains, overall_best) = find_parameters_mcmc_parallel(
        data,
        y_error_bar,
        initial_guess,
        bounds,
        fix_tau,
        num_steps,
        burn_in,
        inwards,
        num_chains,
    );

    // Combine chains for analysis
    let (combined_chain, combined_likelihoods) = combine_chains(&chains);

    // Check convergence
    let r_hat = calculate_gelman_rubin(&chains, inwards);
    println!("Gelman-Rubin R-hat statistics: {:?}", r_hat);
    let converged = r_hat.iter().all(|&r| r < 1.1);

    if !converged {
        println!("Warning: Chains may not have converged (R-hat > 1.1)");
        println!("Consider running more steps or tuning step sizes.");
    } else {
        println!("Chains appear to have converged (R-hat < 1.1)");
    }

    // Create trace plots for each chain
    create_chain_trace_plots(&chains)?;

    Ok((overall_best, combined_chain, combined_likelihoods))
}

#[derive(Clone, Copy)]
struct Walker {
    params: [f64; 4],
    log_prob: f64,
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
) {
    let n = walkers.len();
    let d = 4;

    // Copy walkers so reads are consistent during parallel update
    let walkers_old = walkers.to_vec();

    walkers.par_iter_mut().enumerate().for_each(|(i, walker)| {
        let mut rng = Pcg64::seed_from_u64(1234 + i as u64);

        // choose partner
        let mut j = rng.random_range(0..n);
        while j == i {
            j = rng.random_range(0..n);
        }

        let z = sample_z(&mut rng, a);
        let mut proposal = walker.params;

        for k in 0..d {
            proposal[k] =
                walkers_old[j].params[k] + z * (walker.params[k] - walkers_old[j].params[k]);
        }

        let proposed_log_prob = log_likelihood(proposal);
        let log_accept = (d as f64 - 1.0) * z.ln() + proposed_log_prob - walker.log_prob;

        if log_accept >= 0.0 || rng.random::<f64>() < log_accept.exp() {
            walker.params = proposal;
            walker.log_prob = proposed_log_prob;
        }
    });
}

fn run_ensemble_sampler(
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
    initial_guess: [f64; 4],
    num_steps: usize,
    burn_in: usize,
    n_walkers: usize,
    inwards: bool,
) -> (Vec<[f64; 4]>, Vec<f64>) {
    let mut rng = Pcg64::seed_from_u64(42);
    let a = 2.0;

    // initialize walkers
    let mut walkers = Vec::with_capacity(n_walkers);
    for _ in 0..n_walkers {
        let mut p = initial_guess;
        for i in 0..4 {
            p[i] *= 1.0 + 0.01 * rng.random::<f64>();
        }
        let lp = log_likelihood(p, data, y_error_bar, true, inwards);
        walkers.push(Walker {
            params: p,
            log_prob: lp,
        });
    }

    let mut chain = Vec::new();
    let mut likelihoods = Vec::new();

    for step in 0..num_steps {
        stretch_move_parallel(
            &mut walkers,
            &|p| log_likelihood(p, data, y_error_bar, true, inwards),
            a,
        );

        if step >= burn_in {
            for w in &walkers {
                chain.push(w.params);
                likelihoods.push(w.log_prob);
            }
        }

        if step % 1000 == 0 {
            println!("Step {}", step);
        }
    }

    (chain, likelihoods)
}

pub fn combine_chains(chains: &[([f64; 4], Vec<[f64; 4]>, Vec<f64>)]) -> (Vec<[f64; 4]>, Vec<f64>) {
    let mut combined_chain = Vec::new();
    let mut combined_likelihoods = Vec::new();

    for (_, chain, likelihoods) in chains {
        combined_chain.extend_from_slice(&chain);
        combined_likelihoods.extend_from_slice(&likelihoods);
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

fn find_parameters_mcmc_parallel(
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
    initial_guess: [f64; 4],
    bounds: [(f64, f64); 4],
    fix_tau: Option<f64>,
    num_steps: usize,
    burn_in: usize,
    inwards: bool,
    num_chains: usize,
) -> (Vec<([f64; 4], Vec<[f64; 4]>, Vec<f64>)>, [f64; 4]) {
    println!("Running {} parallel MCMC chains...", num_chains);

    // Create a counter for progress reporting (optional)
    let progress_counter = AtomicUsize::new(0);

    let num_params = match inwards {
        true => 3,
        false => 4,
    };

    // Run chains in parallel
    let chains: Vec<([f64; 4], Vec<[f64; 4]>, Vec<f64>)> = (0..num_chains)
        .into_par_iter()
        .map(|chain_id| {
            // Create a unique seed for each chain
            let seed = 42 + chain_id as u64 * 12345;
            let mut rng = Pcg64::seed_from_u64(seed);

            // Perturb initial guess slightly for each chain
            let mut chain_initial_guess = initial_guess;
            for i in 0..num_params {
                if fix_tau.is_none() || i != 2 {
                    // Add 10% random perturbation
                    chain_initial_guess[i] *= 1.0 + 0.1 * (rng.random::<f64>() - 0.5);
                }
            }

            // Run single chain
            let (best_params, chain, likelihoods) = run_single_chain(
                data,
                y_error_bar,
                chain_initial_guess,
                bounds,
                fix_tau,
                num_steps,
                burn_in,
                inwards,
                &mut rng,
                chain_id,
            );

            // Update progress
            let completed = progress_counter.fetch_add(1, Ordering::Relaxed) + 1;
            println!(
                "Chain {} completed. Progress: {}/{}",
                chain_id, completed, num_chains
            );

            (best_params, chain, likelihoods)
        })
        .collect();

    // Find overall best parameters (from all chains)
    let overall_best = find_overall_best(&chains);

    (chains, overall_best)
}

fn find_overall_best(chains: &[([f64; 4], Vec<[f64; 4]>, Vec<f64>)]) -> [f64; 4] {
    let mut best_log_likelihood = f64::NEG_INFINITY;
    let mut best_params = [0.0; 4];

    for (params, _, likelihoods) in chains {
        if let Some(max_likelihood) = likelihoods.iter().max_by(|a, b| a.partial_cmp(b).unwrap()) {
            let log_likelihood = max_likelihood.ln();
            if log_likelihood > best_log_likelihood {
                best_log_likelihood = log_likelihood;
                best_params = *params;
            }
        }
    }

    best_params
}

pub fn calculate_gelman_rubin(
    chains: &[([f64; 4], Vec<[f64; 4]>, Vec<f64>)],
    inwards: bool,
) -> Vec<f64> {
    let m = chains.len() as f64;
    let n = chains[0].1.len() as f64; // Samples per chain after burn-in

    let num_params = match inwards {
        true => 3,
        false => 4,
    };
    let mut r_hat = vec![0.0; num_params];

    for param_idx in 0..num_params {
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

        dbg!(param_idx);
        dbg!(overall_mean);
        dbg!(within_var);
        dbg!(between_var);
        // Calculate pooled variance
        let pooled_var = ((n - 1.0) / n) * within_var + between_var;

        // R-hat statistic (should approach 1.0 as chains converge)
        r_hat[param_idx] = (pooled_var / within_var); //.sqrt();
    }

    r_hat
}

const LOG_SCALE: [bool; 4] = [true, false, false, true];
const LOG_PRIOR: [bool; 4] = [true, false, false, true];
#[derive(Clone, Copy)]
struct Proposal {
    params: [f64; 4],
    log_jacobian: f64,
    valid: bool,
}

fn propose(
    current: &[f64; 4],
    step_size: &[f64; 4],
    bounds: &[(f64, f64); 4],
    log_scale: &[bool; 4],
    fix_tau: Option<f64>,
    num_params: usize,
    rng: &mut impl Rng,
) -> Proposal {
    let mut params = *current;
    let mut log_jacobian = 0.0;

    for i in 0..num_params {
        // Handle fixed tau
        if fix_tau.is_some() && i == 2 {
            params[i] = fix_tau.unwrap();
            continue;
        }

        let (cur_val, is_log) = if log_scale[i] {
            (current[i].ln(), true)
        } else {
            (current[i], false)
        };

        let normal = Normal::new(cur_val, step_size[i]).unwrap();
        let mut prop = normal.sample(rng);

        if is_log {
            prop = prop.exp();
        }

        // Jacobian to map log sampling back to linear prior
        if !LOG_PRIOR[i] {
            log_jacobian += prop.ln() - current[i].ln();
        }

        // Hard reject if OOB
        if prop < bounds[i].0 || prop > bounds[i].1 {
            return Proposal {
                params: *current,
                log_jacobian: 0.0,
                valid: false,
            };
        }

        params[i] = prop;
    }

    Proposal {
        params,
        log_jacobian,
        valid: true,
    }
}

fn run_single_chain(
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
    initial_guess: [f64; 4], // m200, c200, tau rho_c
    bounds: [(f64, f64); 4],
    fix_tau: Option<f64>,
    num_steps: usize,
    burn_in: usize,
    inwards: bool,
    rng: &mut impl Rng,
    chain_id: usize,
) -> ([f64; 4], Vec<[f64; 4]>, Vec<f64>) {
    let mut current_params = initial_guess;
    let mut current_log_likelihood =
        log_likelihood(current_params, data, y_error_bar, true, inwards);

    let mut accepted = 0;
    let mut chain = Vec::with_capacity(num_steps);
    let mut likelihoods = Vec::with_capacity(num_steps);

    let mut acceptance_history = vec![false; burn_in];

    let num_params = match inwards {
        true => 3,
        false => 4,
    };

    let mut step_size: [f64; 4] = [0.0; 4];
    for i in 0..4 {
        if LOG_SCALE[i] {
            step_size[i] = 0.01 * initial_guess[i].ln()
        } else {
            step_size[i] = 0.01
        }
    }

    for step in 0..num_steps {
        // Propose new parameters
        let proposal = propose(
            &current_params,
            &step_size,
            &bounds,
            &LOG_SCALE,
            fix_tau,
            num_params,
            rng,
        );

        // Calculate log likelihood for proposed parameters
        let proposed_log_likelihood =
            log_likelihood(proposal.params, data, y_error_bar, true, inwards);

        // Metropolis-Hast acceptance ratio
        let log_acceptance =
            proposed_log_likelihood - current_log_likelihood + proposal.log_jacobian;

        if proposal.valid && (log_acceptance >= 0.0 || rng.random::<f64>() < log_acceptance.exp()) {
            current_params = proposal.params;
            current_log_likelihood = proposed_log_likelihood;
            accepted += 1;

            if step < burn_in {
                acceptance_history[step] = true;
            }
        }
        // Adjust step
        const TARGET_RATIO: f64 = 0.234; //Roberts, G. O., Gelman, A., & Gilks, W. R. (1997)
        if step % 10 == 0 && step < burn_in {
            let acceptance_rate = acceptance_history[step.saturating_sub(100)..]
                .iter()
                .filter(|&&b| b)
                .count() as f64
                / (step - step.saturating_sub(100)) as f64;
            // Target acceptance rate between 0.4 and 0.6
            for i in 0..num_params {
                if acceptance_rate < TARGET_RATIO {
                    step_size[i] *= 0.99
                } else if acceptance_rate > TARGET_RATIO {
                    step_size[i] *= 1.01
                }
            }
        }

        // Store after burn-in period
        if step >= burn_in {
            chain.push(current_params);
            likelihoods.push(current_log_likelihood.exp()); // Convert back from log
        }

        // Optional: Print progress
        if step % 1000 == 0 {
            println!(
                "Chain {}: Step {}/{}: Accepted rate: {:.2}%, LogLik: {:.4} \nStep sizes: {}, {}, {}, {}",
                chain_id,
                step,
                num_steps,
                (accepted as f64 / (step + 1) as f64) * 100.0,
                current_log_likelihood,
                step_size[0],
                step_size[1],
                step_size[2],
                step_size[3],
            );
        }
    }

    // Find best parameters (maximum likelihood)
    let best_idx = likelihoods
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(idx, _)| idx)
        .unwrap_or(0);

    let best_params = chain[best_idx];

    println!(
        "Chain {} finished: Acceptance: {:.1}%, Best log-likelyhood: {:.4}",
        chain_id,
        (accepted as f64 / num_steps as f64) * 100.0,
        likelihoods[best_idx].ln()
    );

    (best_params, chain, likelihoods)
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

    let temp = UVB_TEMP;

    while n < 10 || (error[error.len() - 1] > 0.01 && rel_error_change > 0.0 && n < 2000) {
        error.push(get_rms_err_of_fit(temp, params, &data, false, inwards));
        rel_error_change =
            ((error[error.len() - 1] - error[error.len() - 2]) / error[error.len() - 2]).abs();

        let mut grad: [f64; 4] = [0.0; 4];

        let grad_scale = (0.1 * (1.0_f64 - 5e-3).powi(n)).max(1e-3);
        for i in 0..params.len() {
            let mut increased_params = params.clone();
            increased_params[i] += 0.5 * grad_scale * params[i];
            let increase_error = get_rms_err_of_fit(temp, increased_params, &data, false, inwards);

            let mut decreased_params = params.clone();
            decreased_params[i] -= 0.5 * grad_scale * params[i];
            let decrease_error = get_rms_err_of_fit(temp, decreased_params, &data, false, inwards);
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

    plot_function(
        &(0..n).map(|n| n as f64).collect(),
        &error,
        "fit_error.png",
        "Fitting Error",
        "Steps",
        "Error",
        None,
        None,
    )
    .unwrap();

    params
}

fn log_likelihood(
    mut params: [f64; 4],
    data: &Vec<(f64, f64)>,
    y_error_bar: &Vec<(f64, f64)>,
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
        UVB_TEMP,
        params[0],
        params[1],
        params[2],
        rho_c,
        (0.1 * data[0].0, halo.r_crit()),
        false,
    );

    let mut chi_squared = 0.0;
    let mut sum_ln_sigma = 0.0;
    for i in 0..data.len() {
        let point = data[i];

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

        let diff = fit.1[low] - point.1;
        let spread = match diff > 0.0 {
            true => y_error_bar[i].1 - point.1,
            false => point.1 - y_error_bar[i].0,
        };
        chi_squared += (diff / spread).powi(2);
        sum_ln_sigma += spread.ln();
    }

    -0.5 * chi_squared - sum_ln_sigma
}

fn get_rms_err_of_fit(
    temp: f64,
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
        temp,
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
