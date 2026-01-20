use crate::{hydrostatics::isothermal_core_collapse_background, plotting::plot_function};

pub fn find_parameters(data: &Vec<(f64, f64)>, initial_guess: [f64; 4], fix_tau: Option<f64>) -> [f64; 4] {
    let mut params = initial_guess;

    let mut error = Vec::new();
    error.push(1.0);
    error.push(1.0);
    let mut rel_error_change = f64::MAX;
    let mut n = 0;

    let temp = 1e4;

    while n < 10 || (error[error.len() - 1] > 0.01 && rel_error_change > 0.0 && n < 2000) {

        error.push(get_err_of_fit(temp, params, &data));
        rel_error_change = ((error[error.len() - 1] - error[error.len() - 2])/error[error.len() - 2]).abs();

        let mut grad: [f64; 4] = [0.0; 4];

        let grad_scale = (0.1 * (1.0_f64 - 5e-3).powi(n)).max(1e-3);
        for i in 0..params.len() {
            let mut increased_params = params.clone();
            increased_params[i] += 0.5*grad_scale*params[i];
            let increase_error = get_err_of_fit(temp, increased_params, &data);
            
            let mut decreased_params = params.clone();
            decreased_params[i] -= 0.5*grad_scale*params[i];
            let decrease_error = get_err_of_fit(temp, decreased_params, &data);

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

fn get_err_of_fit(temp: f64, params: [f64; 4], data: &Vec<(f64, f64)>) -> f64 {
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