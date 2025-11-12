"""
Random Dot Product Graph (RDPG) Simulation with Spectral Embedding

This script:
1. Generates a network under the Random Dot Product Graph (RDPG) model
2. Uses spectral embedding to estimate the probability matrix P_ij
3. Calculates the Mean Squared Error (MSE) between estimated and true P_ij
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import svd
from typing import Tuple


def generate_latent_positions(n: int, d: int, distribution: str = 'uniform') -> np.ndarray:
    """
    Generate latent positions for nodes.

    Parameters:
    -----------
    n : int
        Number of nodes
    d : int
        Dimension of latent space
    distribution : str
        Distribution type ('uniform' or 'normal')

    Returns:
    --------
    X : np.ndarray
        Latent position matrix of shape (n, d)
    """
    if distribution == 'uniform':
        # Generate positions uniformly in [0, 1]^d
        X = np.random.uniform(0, 1, size=(n, d))
    elif distribution == 'normal':
        # Generate positions from standard normal and normalize
        X = np.random.randn(n, d)
        X = np.abs(X)  # Make positive for valid probabilities
        # Normalize to ensure P_ij <= 1
        X = X / (np.sqrt(d) * 3)  # Scale down to avoid probabilities > 1
    else:
        raise ValueError(f"Unknown distribution: {distribution}")

    return X


def compute_probability_matrix(X: np.ndarray) -> np.ndarray:
    """
    Compute the true probability matrix P from latent positions.

    P_ij = X_i^T * X_j

    Parameters:
    -----------
    X : np.ndarray
        Latent position matrix of shape (n, d)

    Returns:
    --------
    P : np.ndarray
        Probability matrix of shape (n, n)
    """
    P = X @ X.T
    # Clip probabilities to [0, 1] range
    P = np.clip(P, 0, 1)
    # Set diagonal to 0 (no self-loops)
    np.fill_diagonal(P, 0)

    return P


def generate_rdpg_network(P: np.ndarray, symmetric: bool = True) -> np.ndarray:
    """
    Generate adjacency matrix under RDPG model.

    Parameters:
    -----------
    P : np.ndarray
        Probability matrix of shape (n, n)
    symmetric : bool
        If True, generate undirected graph

    Returns:
    --------
    A : np.ndarray
        Adjacency matrix of shape (n, n)
    """
    n = P.shape[0]
    A = np.random.binomial(1, P)

    if symmetric:
        # Make undirected by symmetrizing
        A = np.triu(A, k=1)  # Upper triangular
        A = A + A.T  # Make symmetric

    # No self-loops
    np.fill_diagonal(A, 0)

    return A


def spectral_embedding(A: np.ndarray, d: int, method: str = 'ASE') -> np.ndarray:
    """
    Perform spectral embedding on adjacency matrix.

    Parameters:
    -----------
    A : np.ndarray
        Adjacency matrix of shape (n, n)
    d : int
        Embedding dimension
    method : str
        'ASE' for Adjacency Spectral Embedding

    Returns:
    --------
    X_hat : np.ndarray
        Estimated latent positions of shape (n, d)
    """
    n = A.shape[0]

    if method == 'ASE':
        # Adjacency Spectral Embedding
        # Perform SVD on adjacency matrix
        U, S, Vt = svd(A, full_matrices=False)

        # Take top d singular values/vectors
        U_d = U[:, :d]
        S_d = np.diag(S[:d])

        # Embed: X_hat = U_d * sqrt(S_d)
        X_hat = U_d @ np.sqrt(S_d)

    else:
        raise ValueError(f"Unknown method: {method}")

    return X_hat


def estimate_probability_matrix(X_hat: np.ndarray) -> np.ndarray:
    """
    Estimate probability matrix from embedded positions.

    P_hat_ij = X_hat_i^T * X_hat_j

    Parameters:
    -----------
    X_hat : np.ndarray
        Estimated latent positions of shape (n, d)

    Returns:
    --------
    P_hat : np.ndarray
        Estimated probability matrix of shape (n, n)
    """
    P_hat = X_hat @ X_hat.T
    # Clip to [0, 1]
    P_hat = np.clip(P_hat, 0, 1)
    # Set diagonal to 0
    np.fill_diagonal(P_hat, 0)

    return P_hat


def calculate_mse(P_true: np.ndarray, P_hat: np.ndarray) -> float:
    """
    Calculate Mean Squared Error between true and estimated probability matrices.

    Parameters:
    -----------
    P_true : np.ndarray
        True probability matrix
    P_hat : np.ndarray
        Estimated probability matrix

    Returns:
    --------
    mse : float
        Mean squared error
    """
    # Calculate MSE only for off-diagonal elements (excluding self-loops)
    n = P_true.shape[0]
    mask = ~np.eye(n, dtype=bool)

    mse = np.mean((P_true[mask] - P_hat[mask])**2)

    return mse


def run_simulation(n: int = 100, d: int = 2, distribution: str = 'uniform',
                   seed: int = None) -> dict:
    """
    Run complete RDPG simulation.

    Parameters:
    -----------
    n : int
        Number of nodes
    d : int
        Latent dimension
    distribution : str
        Distribution for latent positions
    seed : int
        Random seed for reproducibility

    Returns:
    --------
    results : dict
        Dictionary containing all results
    """
    if seed is not None:
        np.random.seed(seed)

    print(f"Running RDPG Simulation:")
    print(f"  Number of nodes (n): {n}")
    print(f"  Latent dimension (d): {d}")
    print(f"  Distribution: {distribution}")
    print("-" * 50)

    # Step 1: Generate latent positions
    print("Step 1: Generating latent positions...")
    X = generate_latent_positions(n, d, distribution)

    # Step 2: Compute true probability matrix
    print("Step 2: Computing true probability matrix P...")
    P_true = compute_probability_matrix(X)

    # Step 3: Generate network (adjacency matrix)
    print("Step 3: Generating network under RDPG model...")
    A = generate_rdpg_network(P_true, symmetric=True)
    edge_count = np.sum(A) / 2  # Divide by 2 for undirected
    density = edge_count / (n * (n - 1) / 2)
    print(f"  Network generated: {int(edge_count)} edges, density = {density:.4f}")

    # Step 4: Spectral embedding
    print("Step 4: Performing spectral embedding...")
    X_hat = spectral_embedding(A, d, method='ASE')

    # Step 5: Estimate probability matrix
    print("Step 5: Estimating probability matrix P_hat...")
    P_hat = estimate_probability_matrix(X_hat)

    # Step 6: Calculate MSE
    print("Step 6: Calculating MSE...")
    mse = calculate_mse(P_true, P_hat)
    print(f"  MSE = {mse:.6f}")

    print("-" * 50)
    print("Simulation complete!")

    results = {
        'X': X,
        'P_true': P_true,
        'A': A,
        'X_hat': X_hat,
        'P_hat': P_hat,
        'mse': mse,
        'edge_count': edge_count,
        'density': density
    }

    return results


def plot_results(results: dict, save_path: str = None):
    """
    Visualize simulation results.

    Parameters:
    -----------
    results : dict
        Results from run_simulation
    save_path : str
        Path to save the figure (optional)
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Plot 1: True latent positions (if d=2)
    if results['X'].shape[1] == 2:
        axes[0, 0].scatter(results['X'][:, 0], results['X'][:, 1], alpha=0.6)
        axes[0, 0].set_title('True Latent Positions')
        axes[0, 0].set_xlabel('Dimension 1')
        axes[0, 0].set_ylabel('Dimension 2')
    else:
        axes[0, 0].text(0.5, 0.5, f'Latent dim = {results["X"].shape[1]}',
                        ha='center', va='center')
        axes[0, 0].set_title('True Latent Positions')

    # Plot 2: Adjacency matrix
    im1 = axes[0, 1].imshow(results['A'], cmap='binary', interpolation='nearest')
    axes[0, 1].set_title(f'Adjacency Matrix\n({int(results["edge_count"])} edges)')
    plt.colorbar(im1, ax=axes[0, 1])

    # Plot 3: True probability matrix
    im2 = axes[0, 2].imshow(results['P_true'], cmap='viridis', vmin=0, vmax=1)
    axes[0, 2].set_title('True Probability Matrix P')
    plt.colorbar(im2, ax=axes[0, 2])

    # Plot 4: Estimated latent positions (if d=2)
    if results['X_hat'].shape[1] == 2:
        axes[1, 0].scatter(results['X_hat'][:, 0], results['X_hat'][:, 1], alpha=0.6)
        axes[1, 0].set_title('Estimated Latent Positions')
        axes[1, 0].set_xlabel('Dimension 1')
        axes[1, 0].set_ylabel('Dimension 2')
    else:
        axes[1, 0].text(0.5, 0.5, f'Latent dim = {results["X_hat"].shape[1]}',
                        ha='center', va='center')
        axes[1, 0].set_title('Estimated Latent Positions')

    # Plot 5: Estimated probability matrix
    im3 = axes[1, 1].imshow(results['P_hat'], cmap='viridis', vmin=0, vmax=1)
    axes[1, 1].set_title('Estimated Probability Matrix P_hat')
    plt.colorbar(im3, ax=axes[1, 1])

    # Plot 6: Scatter plot of P_true vs P_hat
    n = results['P_true'].shape[0]
    mask = ~np.eye(n, dtype=bool)
    axes[1, 2].scatter(results['P_true'][mask], results['P_hat'][mask],
                       alpha=0.3, s=1)
    axes[1, 2].plot([0, 1], [0, 1], 'r--', label='Perfect estimation')
    axes[1, 2].set_xlabel('True P_ij')
    axes[1, 2].set_ylabel('Estimated P_ij')
    axes[1, 2].set_title(f'MSE = {results["mse"]:.6f}')
    axes[1, 2].legend()
    axes[1, 2].grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    plt.show()


def main():
    """
    Main function to run the simulation.
    """
    # Run simulation with default parameters
    results = run_simulation(n=100, d=2, distribution='uniform', seed=42)

    # Plot results
    plot_results(results, save_path='rdpg_simulation_results.png')

    # Run multiple simulations to see MSE variation
    print("\n" + "="*50)
    print("Running multiple simulations to analyze MSE...")
    print("="*50 + "\n")

    n_simulations = 10
    mse_values = []

    for i in range(n_simulations):
        results_i = run_simulation(n=100, d=2, distribution='uniform', seed=i)
        mse_values.append(results_i['mse'])
        print(f"Simulation {i+1}/{n_simulations}: MSE = {results_i['mse']:.6f}")

    print("\n" + "-"*50)
    print(f"Average MSE across {n_simulations} simulations: {np.mean(mse_values):.6f}")
    print(f"Std Dev of MSE: {np.std(mse_values):.6f}")
    print(f"Min MSE: {np.min(mse_values):.6f}")
    print(f"Max MSE: {np.max(mse_values):.6f}")


if __name__ == "__main__":
    main()
