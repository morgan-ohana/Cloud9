
use rand_distr::{Distribution, Normal};
use rand::Rng;

use crate::{hydrostatics::isothermal_core_collapse_background, plotting::plot_function};

pub fn find_parameters_mcmc(
    data: &Vec<(f64, f64)>,
    initial_guess: [f64; 4],
    fix_tau: Option<f64>,
    num_steps: usize,
    burn_in: usize,
) -> ([f64; 4], Vec<[f64; 4]>, Vec<f64>) {
    let mut current_params = initial_guess;
    let mut current_log_likelihood = log_likelihood(current_params, data);
    
    let mut accepted = 0;
    let mut chain = Vec::with_capacity(num_steps);
    let mut likelihoods = Vec::with_capacity(num_steps);
    
    // Normal distributions for proposals (centered at current values)
    let mut rng = rand::rng();
    
    let mut relative_step_size = 0.01;
    for step in 0..num_steps {
        // Propose new parameters
        let mut proposed_params = current_params;
        for i in 0..4 {
            // Skip fixed parameter (tau is index 2)
            if fix_tau.is_some() && i == 2 {
                proposed_params[i] = fix_tau.unwrap();
                continue;
            }
            
            let normal = Normal::new(current_params[i], relative_step_size * current_params[i].abs()).unwrap();
            proposed_params[i] = normal.sample(&mut rng);
            
            // Ensure parameters stay positive if needed
            if proposed_params[i] < 0.0 {
                proposed_params[i] = 1e-10;
            }
        }
        
        // Calculate log likelihood for proposed parameters
        let proposed_log_likelihood = log_likelihood(proposed_params, data);
        
        // Metropolis-Hast acceptance ratio
        let acceptance_ratio = (proposed_log_likelihood - current_log_likelihood).exp();
        
        if acceptance_ratio >= 1.0 || rng.random::<f64>() < acceptance_ratio {
            current_params = proposed_params;
            current_log_likelihood = proposed_log_likelihood;
            accepted += 1;
        }

        // Adjust step
        if step % 100 == 0 {
            let acceptance_rate = accepted as f64 / (step + 1) as f64;
            // Target acceptance rate between 0.4 and 0.6
            if acceptance_rate < 0.4 {
                relative_step_size *= 0.99;
            } else if acceptance_rate > 0.6 {
                relative_step_size *= 1.01;
            }
        }

        // Store after burn-in period
        if step >= burn_in {
            chain.push(current_params);
            likelihoods.push(current_log_likelihood.exp()); // Convert back from log
        }
        
        // Optional: Print progress
        if step % 1000 == 0 {
            println!("Step {}/{}: Accepted rate: {:.2}%, LogLik: {:.4}, RelStep: {:.4}", 
                step, num_steps, 
                (accepted as f64 / (step + 1) as f64) * 100.0,
                current_log_likelihood,
                relative_step_size
            );
        }
    }
    
    // Calculate final acceptance rate
    let acceptance_rate = accepted as f64 / num_steps as f64;
    println!("Final acceptance rate: {:.2}%", acceptance_rate * 100.0);
    
    // Find best parameters (maximum likelihood)
    let best_idx = likelihoods
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(idx, _)| idx)
        .unwrap_or(0);
    
    let best_params = chain[best_idx];
    
    // Optional: Calculate parameter statistics
    calculate_statistics(&chain, &best_params);
    
    (best_params, chain, likelihoods)
}

pub fn find_parameters_gradient_descent(data: &Vec<(f64, f64)>, initial_guess: [f64; 4], fix_tau: Option<f64>) -> [f64; 4] {
    let mut params = initial_guess;

    let mut error = Vec::new();
    error.push(1.0);
    error.push(1.0);
    let mut rel_error_change = f64::MAX;
    let mut n = 0;

    let temp = 1e4;

    while n < 10 || (error[error.len() - 1] > 0.01 && rel_error_change > 0.0 && n < 2000) {

        error.push(get_rms_err_of_fit(temp, params, &data));
        rel_error_change = ((error[error.len() - 1] - error[error.len() - 2])/error[error.len() - 2]).abs();

        let mut grad: [f64; 4] = [0.0; 4];

        let grad_scale = (0.1 * (1.0_f64 - 5e-3).powi(n)).max(1e-3);
        for i in 0..params.len() {
            let mut increased_params = params.clone();
            increased_params[i] += 0.5*grad_scale*params[i];
            let increase_error = get_rms_err_of_fit(temp, increased_params, &data);
            
            let mut decreased_params = params.clone();
            decreased_params[i] -= 0.5*grad_scale*params[i];
            let decrease_error = get_rms_err_of_fit(temp, decreased_params, &data);
            grad[i] = (increase_error - decrease_error) / (grad_scale*params[i])
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

    plot_function(&(0..n).map(|n| n as f64).collect(),
        &error,
        "fit_error.png",
        "Fitting Error",
        "Steps", 
        "Error",
        None,
        None
    ).unwrap();

    params
}

fn log_likelihood(params: [f64; 4], data: &Vec<(f64, f64)>) -> f64 {
    let error = get_rms_err_of_fit(1e4, params, data);
    
    // Gaussian log-likelihood: -0.5 * χ²
    // where χ² = n * error² (since error is RMS)
    let n = data.len() as f64;
    let chi_squared = n * error.powi(2);
    
    -0.5 * chi_squared
}

fn get_rms_err_of_fit(temp: f64, params: [f64; 4], data: &Vec<(f64, f64)>) -> f64 {
    let fit = isothermal_core_collapse_background(
        temp,
        params[0],
        params[1],
        params[2],
        (0.1*data[0].0, 10.0*data.last().unwrap().0),
        params[3],
        false
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

        error += ((fit.1[low] - point.1)/point.1).powi(2);
    }

    error /= data.len() as f64;
    error.sqrt()
}

fn calculate_statistics(chain: &[[f64; 4]], best_params: &[f64; 4]) {
    let n = chain.len() as f64;
    
    println!("\nParameter Statistics:");
    println!("{:<10} {:<12} {:<12} {:<12}", "Param", "Best", "Mean", "Std Dev");
    println!("{}", "-".repeat(48));

    let names = ["rho_s_0", "r_s_0", "tau", "rho_c"];
    
    for i in 0..4 {
        let mean = chain.iter().map(|p| p[i]).sum::<f64>() / n;
        let variance = chain.iter()
            .map(|p| (p[i] - mean).powi(2))
            .sum::<f64>() / n;
        let std_dev = variance.sqrt();
        
        println!("{:<10} {:<12} {:<12} {:<12}", 
                 names[i], 
                 format_number(best_params[i]), 
                 format_number(mean), 
                 format_number(std_dev));
    }
}

fn format_number(x: f64) -> String {
    // Use scientific notation for very large or very small numbers
    if x.abs() >= 1e3 || (x.abs() > 0.0 && x.abs() < 1e-2) {
        format!("{:.4e}", x)
    } else {
        format!("{:.4}", x)
    }
}