# NCBI - Random Dot Product Graph Simulation

This repository contains a simulation for Random Dot Product Graph (RDPG) model with spectral embedding for probability matrix estimation.

## Overview

The simulation performs the following steps:

1. **Generate Network Data**: Creates a network under the Random Dot Product Graph (RDPG) model
2. **Spectral Embedding**: Uses Adjacency Spectral Embedding (ASE) to estimate latent positions
3. **Probability Matrix Estimation**: Estimates P_ij from the embedded positions
4. **MSE Calculation**: Computes Mean Squared Error between true and estimated probability matrices

## Random Dot Product Graph Model

In the RDPG model:
- Each node i has a latent position vector X_i ∈ ℝ^d
- The probability of an edge between nodes i and j is: **P_ij = X_i^T · X_j**
- An adjacency matrix A is generated where A_ij ~ Bernoulli(P_ij)

## Spectral Embedding

The Adjacency Spectral Embedding (ASE) method:
1. Performs SVD on the adjacency matrix: A = UΣV^T
2. Estimates latent positions: X̂ = U_d · √(Σ_d), where d is the embedding dimension
3. Estimates probability matrix: P̂_ij = X̂_i^T · X̂_j

## Installation

Install the required dependencies:

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

Run the simulation with default parameters:

```bash
python rdpg_simulation.py
```

This will:
- Generate a network with 100 nodes and 2-dimensional latent space
- Perform spectral embedding
- Calculate MSE
- Generate visualization plots
- Run 10 simulations to analyze MSE variation

### Custom Parameters

You can modify the simulation parameters in the `main()` function or use the `run_simulation()` function directly:

```python
from rdpg_simulation import run_simulation, plot_results

# Run simulation with custom parameters
results = run_simulation(
    n=200,           # Number of nodes
    d=3,             # Latent dimension
    distribution='uniform',  # 'uniform' or 'normal'
    seed=42          # Random seed
)

# Plot results
plot_results(results)
```

## Output

The simulation produces:

1. **Console Output**:
   - Number of edges and network density
   - MSE between true and estimated probability matrices
   - Statistics across multiple runs

2. **Visualization** (`rdpg_simulation_results.png`):
   - True latent positions (if d=2)
   - Adjacency matrix
   - True probability matrix P
   - Estimated latent positions (if d=2)
   - Estimated probability matrix P̂
   - Scatter plot comparing P_ij and P̂_ij

## Key Functions

- `generate_latent_positions()`: Generate latent positions for nodes
- `compute_probability_matrix()`: Compute true P from latent positions
- `generate_rdpg_network()`: Generate adjacency matrix under RDPG
- `spectral_embedding()`: Perform ASE on adjacency matrix
- `estimate_probability_matrix()`: Estimate P̂ from embedded positions
- `calculate_mse()`: Calculate MSE between true and estimated P
- `run_simulation()`: Run complete simulation pipeline
- `plot_results()`: Visualize simulation results

## References

- Sussman, D. L., Tang, M., Fishkind, D. E., & Priebe, C. E. (2012). A consistent adjacency spectral embedding for stochastic blockmodel graphs. Journal of the American Statistical Association, 107(499), 1119-1128.
- Athreya, A., Priebe, C. E., Tang, M., Lyzinski, V., Marchette, D. J., & Sussman, D. L. (2016). A limit theorem for scaled eigenvectors of random dot product graphs. Sankhya A, 78(1), 1-18.