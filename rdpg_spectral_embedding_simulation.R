# ============================================================================
# Random Dot Product Graph (RDPG) with Spectral Embedding Simulation
# ============================================================================
# This script:
# 1. Generates a network under the Random Dot Product Graph model
# 2. Uses spectral embedding to estimate P_ij (edge probabilities)
# 3. Calculates MSE between estimated and true P_ij

# Load required libraries
library(Matrix)
library(irlba)  # for efficient eigendecomposition

# Set seed for reproducibility
set.seed(123)

# ============================================================================
# PARAMETERS
# ============================================================================
n <- 100          # Number of nodes
d <- 3            # Dimension of latent space
eigen_dim <- 3    # Dimension for spectral embedding (should equal d)

# ============================================================================
# STEP 1: Generate Random Dot Product Graph
# ============================================================================

cat("Step 1: Generating RDPG network...\n")

# Generate latent positions X (n x d matrix)
# Each row is a latent position for a node
# We'll use uniform distribution in [0,1]^d to ensure valid probabilities
X <- matrix(runif(n * d, 0, 0.5), nrow = n, ncol = d)

# Normalize to ensure P_ij values are valid probabilities
# (optional, but helps ensure X_i^T X_j <= 1)
X <- t(apply(X, 1, function(x) x / sqrt(sum(x^2)) * 0.7))

# Calculate true probability matrix P
# P_ij = X_i^T X_j (inner product of latent positions)
P_true <- X %*% t(X)

# Ensure P is symmetric and bounded [0,1]
P_true <- pmin(pmax(P_true, 0), 1)
diag(P_true) <- 0  # No self-loops

cat(sprintf("  - Generated %d x %d probability matrix\n", n, n))
cat(sprintf("  - P_ij range: [%.4f, %.4f]\n", min(P_true), max(P_true)))

# Sample adjacency matrix A from P_true
# A_ij ~ Bernoulli(P_ij)
A <- matrix(0, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if (runif(1) < P_true[i, j]) {
      A[i, j] <- 1
      A[j, i] <- 1  # Make symmetric
    }
  }
}

num_edges <- sum(A) / 2
cat(sprintf("  - Generated network with %d edges (density: %.4f)\n",
            num_edges, num_edges / (n * (n-1) / 2)))

# ============================================================================
# STEP 2: Spectral Embedding
# ============================================================================

cat("\nStep 2: Performing spectral embedding...\n")

# Compute eigendecomposition of adjacency matrix
# For RDPG, we use the top d eigenvectors
eigen_result <- eigen(A, symmetric = TRUE)

# Extract top d eigenvectors and eigenvalues
eigenvalues <- eigen_result$values[1:eigen_dim]
eigenvectors <- eigen_result$vectors[, 1:eigen_dim]

cat(sprintf("  - Top %d eigenvalues: %s\n", eigen_dim,
            paste(round(eigenvalues, 3), collapse=", ")))

# Estimate latent positions: X_hat = U * sqrt(Lambda)
# where U is the matrix of eigenvectors and Lambda is diagonal matrix of eigenvalues
Lambda_sqrt <- diag(sqrt(pmax(eigenvalues, 0)))  # Take positive part
X_hat <- eigenvectors %*% Lambda_sqrt

cat(sprintf("  - Estimated latent positions: %d x %d matrix\n", nrow(X_hat), ncol(X_hat)))

# ============================================================================
# STEP 3: Estimate P_ij from Spectral Embedding
# ============================================================================

cat("\nStep 3: Estimating P_ij from spectral embedding...\n")

# Calculate estimated probability matrix
# P_hat_ij = X_hat_i^T X_hat_j
P_hat <- X_hat %*% t(X_hat)

# Ensure P_hat is bounded [0,1] and symmetric
P_hat <- pmin(pmax(P_hat, 0), 1)
diag(P_hat) <- 0  # No self-loops

cat(sprintf("  - P_hat range: [%.4f, %.4f]\n", min(P_hat), max(P_hat)))

# ============================================================================
# STEP 4: Calculate MSE
# ============================================================================

cat("\nStep 4: Calculating MSE between P_true and P_hat...\n")

# Calculate MSE for upper triangular part (to avoid double counting)
upper_tri_idx <- upper.tri(P_true)
P_true_vec <- P_true[upper_tri_idx]
P_hat_vec <- P_hat[upper_tri_idx]

MSE <- mean((P_true_vec - P_hat_vec)^2)
RMSE <- sqrt(MSE)
MAE <- mean(abs(P_true_vec - P_hat_vec))

cat(sprintf("  - MSE: %.6f\n", MSE))
cat(sprintf("  - RMSE: %.6f\n", RMSE))
cat(sprintf("  - MAE: %.6f\n", MAE))

# Calculate correlation between true and estimated probabilities
cor_P <- cor(P_true_vec, P_hat_vec)
cat(sprintf("  - Correlation: %.4f\n", cor_P))

# ============================================================================
# STEP 5: Visualization (optional)
# ============================================================================

cat("\nStep 5: Creating diagnostic plots...\n")

# Create plots if possible
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Plot 1: True vs Estimated P_ij
  plot_data <- data.frame(
    P_true = P_true_vec,
    P_hat = P_hat_vec
  )

  p1 <- ggplot(plot_data, aes(x = P_true, y = P_hat)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = "True vs Estimated Edge Probabilities",
         x = "True P_ij",
         y = "Estimated P_ij") +
    theme_minimal() +
    coord_fixed()

  ggsave("rdpg_true_vs_estimated.png", p1, width = 6, height = 6)
  cat("  - Saved plot: rdpg_true_vs_estimated.png\n")

  # Plot 2: Distribution of errors
  plot_data$error <- plot_data$P_hat - plot_data$P_true

  p2 <- ggplot(plot_data, aes(x = error)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Estimation Errors",
         x = "Error (P_hat - P_true)",
         y = "Frequency") +
    theme_minimal()

  ggsave("rdpg_error_distribution.png", p2, width = 6, height = 4)
  cat("  - Saved plot: rdpg_error_distribution.png\n")
} else {
  cat("  - ggplot2 not available. Skipping visualization.\n")
  cat("  - Install with: install.packages('ggplot2')\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n" , rep("=", 70), "\n", sep="")
cat("SIMULATION SUMMARY\n")
cat(rep("=", 70), "\n", sep="")
cat(sprintf("Network size: %d nodes\n", n))
cat(sprintf("Latent dimension: %d\n", d))
cat(sprintf("Number of edges: %d\n", num_edges))
cat(sprintf("Network density: %.4f\n", num_edges / (n * (n-1) / 2)))
cat(sprintf("\nPerformance Metrics:\n"))
cat(sprintf("  MSE:  %.6f\n", MSE))
cat(sprintf("  RMSE: %.6f\n", RMSE))
cat(sprintf("  MAE:  %.6f\n", MAE))
cat(sprintf("  Correlation: %.4f\n", cor_P))
cat(rep("=", 70), "\n", sep="")

# Return results
results <- list(
  X_true = X,
  X_hat = X_hat,
  P_true = P_true,
  P_hat = P_hat,
  A = A,
  MSE = MSE,
  RMSE = RMSE,
  MAE = MAE,
  correlation = cor_P
)

cat("\nResults stored in 'results' list.\n")
cat("Access with: results$MSE, results$P_true, results$P_hat, etc.\n")
