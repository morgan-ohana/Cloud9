use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;
use rayon::prelude::*;
use rkyv::{Archive, Deserialize, Serialize};

#[derive(Archive, Deserialize, Serialize, Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct MCMCOutput {
    pub best_params: Vec<f64>,
    pub n_params: usize,
    pub flattened_chain: Vec<f64>,
    pub log_likelihoods: Vec<f64>,
    pub gelman_rubin: Vec<f64>,
}

impl MCMCOutput {
    pub fn init(
        best_params: Vec<f64>,
        chain: Vec<Vec<f64>>,
        log_likelihoods: Vec<f64>,
        gelman_rubin: Vec<f64>,
    ) -> Self {
        MCMCOutput {
            best_params,
            n_params: chain[0].len(),
            flattened_chain: chain.into_iter().flatten().collect(),
            log_likelihoods,
            gelman_rubin,
        }
    }

    pub fn chain(&self) -> Vec<Vec<f64>> {
        self.flattened_chain
            .chunks(self.n_params)
            .map(|c| c.to_vec())
            .collect()
    }
}

pub trait MCMCCore: Sync {
    fn get_bounds(&self) -> &[[f64; 2]];
    fn get_log_likelihood(&self, params: &[f64]) -> f64;
}

#[derive(Clone)]
struct Walker {
    params: Vec<f64>,
    log_prob: f64,
    rng: Pcg64,
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

fn stretch_move_parallel(walkers: &mut [Walker], core: &impl MCMCCore, a: f64) {
    let n = walkers.len();
    let bounds = core.get_bounds();
    let d = bounds.len();

    // Copy walkers so reads are consistent during parallel update
    let walkers_old = walkers.to_vec();

    walkers.par_iter_mut().enumerate().for_each(|(i, walker)| {
        // choose partner
        let mut j = walker.rng.random_range(0..n);
        while j == i {
            j = walker.rng.random_range(0..n);
        }

        let z = sample_z(&mut walker.rng, a);
        let proposal: Vec<f64> = walker
            .params
            .iter()
            .zip(&walkers_old[j].params)
            .map(|(&xi, &xj)| xj + z * (xi - xj))
            .collect();

        let in_bounds = proposal
            .iter()
            .zip(bounds)
            .all(|(&p, b)| p >= b[0] && p <= b[1]);

        if !in_bounds {
            // if out of bounds just stop so walker stays put ie rejecting move
            return;
        }

        let proposed_log_prob = core.get_log_likelihood(&proposal);

        let log_accept = (d as f64 - 1.0) * z.ln() + proposed_log_prob - walker.log_prob;

        if log_accept >= 0.0 || walker.rng.random::<f64>() < log_accept.exp() {
            walker.params = proposal;
            walker.log_prob = proposed_log_prob;
        }
    });
}

pub fn mcmc(
    core: &impl MCMCCore,
    num_steps: usize,
    burn_in: usize,
    n_walkers: usize,
) -> MCMCOutput {
    let mut rng = Pcg64::seed_from_u64(42);
    let a = 2.0;

    // initialize walkers
    let mut walkers: Vec<Walker> = (0..n_walkers)
        .map(|i| {
            let params: Vec<f64> = core
                .get_bounds()
                .iter()
                .map(|b| b[0] + rng.random::<f64>() * (b[1] - b[0]))
                .collect();
            let log_prob = core.get_log_likelihood(&params);

            Walker {
                params,
                log_prob,
                rng: Pcg64::seed_from_u64(42 + i as u64 * 7919),
            }
        })
        .collect();

    let mut chains: Vec<Vec<Vec<f64>>> = vec![Vec::with_capacity(num_steps); n_walkers];
    let mut log_likelihoods: Vec<Vec<f64>> = vec![Vec::with_capacity(num_steps); n_walkers];

    for step in 0..(num_steps + burn_in) {
        stretch_move_parallel(&mut walkers, core, a);

        if step >= burn_in {
            for w in 0..n_walkers {
                chains[w].push(walkers[w].params.clone());
                log_likelihoods[w].push(walkers[w].log_prob);
            }
        }

        if step % 1000 == 0 {
            if step < burn_in {
                println!("Burn in step {step}")
            } else {
                println!("Step {}", step - burn_in);
            }
        }
    }

    // Split each walker's chain in half, giving 2*n_walkers chains
    let half = num_steps / 2;
    let split: Vec<Vec<Vec<f64>>> = chains
        .iter()
        .flat_map(|walker_chain| {
            // truncate to even length so both halves are equal
            let trimmed = &walker_chain[..half * 2];
            vec![trimmed[..half].to_vec(), trimmed[half..].to_vec()]
        })
        .collect();
    let gelman_rubin = calculate_gelman_rubin(&split);

    let combined_chain: Vec<Vec<f64>> = chains.into_iter().flatten().collect();
    let combined_log_likelihoods: Vec<f64> = log_likelihoods.into_iter().flatten().collect();

    // Best params
    let best_idx = combined_log_likelihoods
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let best_params = combined_chain[best_idx].clone();

    MCMCOutput::init(
        best_params,
        combined_chain,
        combined_log_likelihoods,
        gelman_rubin,
    )
}

pub fn calculate_gelman_rubin(chains: &[Vec<Vec<f64>>]) -> Vec<f64> {
    let m = chains.len() as f64; // number of chains
    let n = chains[0].len() as f64; // Samples per chain
    let n_params = chains[0][0].len();

    let mut r_hat = vec![0.0; n_params];

    for param_idx in 0..n_params {
        // Collect samples for this parameter
        let mut param_samples = Vec::new();
        for chain in chains {
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
        r_hat[param_idx] = if within_var > 0.0 {
            (pooled_var / within_var).sqrt()
        } else {
            f64::NAN
        };
    }

    r_hat
}
