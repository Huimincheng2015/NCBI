# RDPG Spectral Embedding Simulation

This R script simulates a Random Dot Product Graph (RDPG) and estimates edge probabilities using spectral embedding.

## What it does

1. **Generates RDPG Network**: Creates latent positions for n nodes in d-dimensional space and generates a network where P(edge between i,j) = X_i^T X_j
2. **Spectral Embedding**: Uses eigendecomposition of the adjacency matrix to estimate latent positions
3. **Estimates P_ij**: Reconstructs edge probability matrix from estimated latent positions
4. **Calculates MSE**: Computes Mean Squared Error between true and estimated probabilities

## Requirements

```r
# Required packages
install.packages("Matrix")
install.packages("irlba")

# Optional (for visualization)
install.packages("ggplot2")
```

## Usage

### Run the simulation

```r
source("rdpg_spectral_embedding_simulation.R")
```

### Access results

The script stores results in a `results` list:

```r
# Access different components
results$MSE           # Mean Squared Error
results$RMSE          # Root Mean Squared Error
results$MAE           # Mean Absolute Error
results$correlation   # Correlation between true and estimated P_ij
results$P_true        # True probability matrix
results$P_hat         # Estimated probability matrix
results$A             # Adjacency matrix
results$X_true        # True latent positions
results$X_hat         # Estimated latent positions
```

## Parameters

You can modify these parameters in the script:

```r
n <- 100          # Number of nodes
d <- 3            # Dimension of latent space
eigen_dim <- 3    # Dimension for spectral embedding
```

## Output

The script produces:
- Console output with step-by-step progress and summary statistics
- `rdpg_true_vs_estimated.png` - Scatter plot of true vs estimated probabilities
- `rdpg_error_distribution.png` - Histogram of estimation errors

## Method

**Random Dot Product Graph Model:**
- Each node i has a latent position X_i ∈ ℝ^d
- Edge probability: P_ij = X_i^T X_j
- Adjacency matrix: A_ij ~ Bernoulli(P_ij)

**Spectral Embedding:**
- Compute eigendecomposition: A = UΛU^T
- Estimate latent positions: X̂ = U_d √Λ_d (top d eigenvectors/values)
- Estimate probabilities: P̂_ij = X̂_i^T X̂_j

## Example Output

```
SIMULATION SUMMARY
======================================================================
Network size: 100 nodes
Latent dimension: 3
Number of edges: 734
Network density: 0.1482

Performance Metrics:
  MSE:  0.002134
  RMSE: 0.046197
  MAE:  0.031245
  Correlation: 0.8765
======================================================================
```
