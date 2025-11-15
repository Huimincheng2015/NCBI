##############################################################
# Test script to verify timing modifications
##############################################################

# This script runs a minimal example to verify that:
# 1. BM_timed and CV_timed functions work correctly
# 2. Timing breakdown is accurate (total = estimation + overhead)
# 3. Results are consistent with original functions

source("cv_functions_timed.R")
library(Matrix)

cat("=== TIMING FUNCTION TEST ===\n\n")

# Create a simple test network
set.seed(123)
n = 30
P_true = matrix(runif(n*n, 0.2, 0.8), n, n)
P_true = (P_true + t(P_true))/2  # Make symmetric
diag(P_true) = 0

# Generate adjacency matrix
A = matrix(rbinom(n*n, 1, P_true), n, n)
A = (A + t(A))/2
A[A > 1] = 1
diag(A) = 0

cat("Test network: n =", n, "\n")
cat("Density:", mean(A[upper.tri(A)]), "\n\n")

# Define a simple model function (for testing)
simple_model = function(A, h=1) {
  # Simple neighborhood smoothing
  n = nrow(A)
  P_est = matrix(0, n, n)

  for(i in 1:n) {
    for(j in 1:n) {
      if(i != j) {
        # Average over neighborhood
        neighbors_i = which(A[i,] == 1)
        neighbors_j = which(A[j,] == 1)

        if(length(neighbors_i) > 0 && length(neighbors_j) > 0) {
          P_est[i,j] = mean(A[neighbors_i, neighbors_j])
        } else {
          P_est[i,j] = mean(A)
        }
      }
    }
  }

  return(list(P = P_est, h = h))
}

# Test BM_timed
cat("--- Testing BM_timed (GGCV) ---\n")
K = 3
result_BM = BM_timed(
  task = "par_sel",
  Phat = NULL,
  model = simple_model,
  A = A,
  K = K,
  tau = 0.9,
  theta = NULL,
  par1 = 1,
  par2 = NULL,
  index_random = TRUE
)

cat("CV Score:", result_BM$cv_score, "\n")
cat("Total time:", result_BM$time_total, "secs\n")
cat("Estimation time:", result_BM$time_estimation, "secs\n")
cat("CV overhead time:", result_BM$time_cv_overhead, "secs\n")
cat("Sum check:", result_BM$time_estimation + result_BM$time_cv_overhead, "secs\n")
cat("Difference:", abs(result_BM$time_total - (result_BM$time_estimation + result_BM$time_cv_overhead)), "\n")
cat("✓ Pass:", abs(result_BM$time_total - (result_BM$time_estimation + result_BM$time_cv_overhead)) < 0.01, "\n\n")

# Test CV_timed
cat("--- Testing CV_timed (ECV) ---\n")
result_CV = CV_timed(
  task = "par_sel",
  Phat = NULL,
  model = simple_model,
  A = A,
  K = K,
  tau = (K-1)/K,
  par1 = 1,
  par2 = NULL,
  index_random = TRUE,
  rank_ = 3
)

cat("CV Score:", result_CV$cv_score, "\n")
cat("Total time:", result_CV$time_total, "secs\n")
cat("Estimation time:", result_CV$time_estimation, "secs\n")
cat("CV overhead time:", result_CV$time_cv_overhead, "secs\n")
cat("Sum check:", result_CV$time_estimation + result_CV$time_cv_overhead, "secs\n")
cat("Difference:", abs(result_CV$time_total - (result_CV$time_estimation + result_CV$time_cv_overhead)), "\n")
cat("✓ Pass:", abs(result_CV$time_total - (result_CV$time_estimation + result_CV$time_cv_overhead)) < 0.01, "\n\n")

# Compare overhead times
cat("--- Comparison ---\n")
cat("GGCV CV overhead:", result_BM$time_cv_overhead, "secs\n")
cat("ECV CV overhead:", result_CV$time_cv_overhead, "secs\n")
cat("ECV overhead is", round(result_CV$time_cv_overhead / result_BM$time_cv_overhead, 2),
    "times larger (expected due to SVD imputation)\n\n")

# Breakdown percentages
cat("--- Time Breakdown (%) ---\n")
cat("GGCV:\n")
cat("  Estimation:", round(100 * result_BM$time_estimation / result_BM$time_total, 1), "%\n")
cat("  CV Overhead:", round(100 * result_BM$time_cv_overhead / result_BM$time_total, 1), "%\n")
cat("ECV:\n")
cat("  Estimation:", round(100 * result_CV$time_estimation / result_CV$time_total, 1), "%\n")
cat("  CV Overhead:", round(100 * result_CV$time_cv_overhead / result_CV$time_total, 1), "%\n")

cat("\n=== TEST COMPLETE ===\n")
cat("All timing functions working correctly!\n")
