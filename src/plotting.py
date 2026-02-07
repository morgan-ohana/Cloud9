import json
import numpy as np
import matplotlib.pyplot as plt
import corner

def load_mcmc_output(json_file):
    """Load MCMC output from JSON file."""
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    best_params = np.array(data['best_params'])
    chain = np.array(data['chain'])
    likelihoods = np.array(data['likelihoods'])
    
    print(f"Loaded MCMC output from {json_file}")
    print(f"Best parameters: {best_params}")
    print(f"Chain shape: {chain.shape}")
    print(f"Likelihoods shape: {likelihoods.shape}")
    
    return best_params, chain, likelihoods

def plot_corner_plot(chain, parameter_names, save_path, truths=None):
    """
    Create a corner plot from MCMC chains.
    
    Parameters:
    -----------
    chain : numpy array of shape (n_samples, n_parameters)
    parameter_names : list of strings, names of parameters
    truths : list or array, true values (if known) or best-fit values
    save_path : str, path to save the plot (if None, displays the plot)
    """
    
    for n in range(len(chain)):
        chain[n, 0] = np.log10(chain[n, 0])
        chain[n, 3] = np.log10(chain[n, 3])

    if truths is not None:
        truths[0] = np.log10(truths[0])
        truths[3] = np.log10(truths[3])
    
    # Create corner plot
    fig = corner.corner(
        chain,
        labels=parameter_names,
        truths=truths,
        truth_color='red',
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 12},
        label_kwargs={"fontsize": 14},
        hist_kwargs={'linewidth': 2},
        plot_datapoints=False,
        fill_contours=False,
        levels=[0.68, 0.95],  # 1-sigma and 2-sigma contours
        smooth=1.0,
        smooth1d=1.0
    )
    
    # Customize plot
    plt.suptitle('MCMC Corner Plot', fontsize=16, y=1.02)
    
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved corner plot to {save_path}")

def compute_statistics(chain):
    """Compute basic statistics from MCMC chains."""
    n_params = chain.shape[1]
    
    print("\nMCMC Chain Statistics:")
    print("=" * 50)
    
    # For each parameter
    for i in range(n_params):
        param_chain = chain[:, i]
        
        # Compute percentiles (16%, 50%, 84% for 1-sigma credible intervals)
        percentiles = np.percentile(param_chain, [16, 50, 84])
        mean = np.mean(param_chain)
        std = np.std(param_chain)
        
        print(f"Parameter {i+1}:")
        print(f"  Mean ± std: {mean:.4f} ± {std:.4f}")
        print(f"  Median (50%): {percentiles[1]:.4f}")
        print(f"  68% credible interval: [{percentiles[0]:.4f}, {percentiles[2]:.4f}]")
        print(f"  Min, Max: [{np.min(param_chain):.4f}, {np.max(param_chain):.4f}]")
        print()

def plot_likelihood_evolution(likelihoods, save_path):
    """Plot the evolution of likelihood during MCMC."""
    plt.figure(figsize=(12, 4))
    
    # Plot full likelihood evolution
    plt.subplot(1, 2, 1)
    plt.plot(likelihoods, alpha=0.7, linewidth=0.5)
    plt.xlabel('MCMC Step')
    plt.ylabel('Likelihood')
    plt.title('Full Likelihood Evolution')
    plt.grid(True, alpha=0.3)
    
    # Plot running best likelihood
    running_best = np.maximum.accumulate(likelihoods)
    plt.subplot(1, 2, 2)
    plt.plot(running_best, 'r-', linewidth=1)
    plt.xlabel('MCMC Step')
    plt.ylabel('Best Likelihood')
    plt.title('Running Best Likelihood')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved likelihood plot to {save_path}")

def main():
    """Main function to load and analyze MCMC output."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python analyze_mcmc.py <json_file>")
        print("Example: python analyze_mcmc.py mcmc_output.json")
        sys.exit(1)
    
    json_file = sys.argv[1]
    
    # Load the data
    best_params, chain, likelihoods = load_mcmc_output(json_file)
    
    # Define parameter names (customize these based on your model)
    parameter_names = [r'log $M_{200}$', r'$c_{200}$', r'$\tau$', r'log $\rho_c$']
    
    # Compute statistics
    compute_statistics(chain)
    
    # Plot corner plot
    print("\nCreating corner plot...")
    plot_corner_plot(
        chain,
        parameter_names=parameter_names,
        save_path='corner_plot_py.png',
        truths=best_params,  # Use best parameters as reference
    )
    
    # Plot likelihood evolution
    print("\nCreating likelihood plots...")
    plot_likelihood_evolution(likelihoods, save_path='likelihood_evolution.png')
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()