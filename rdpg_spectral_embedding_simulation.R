set.seed(123)

n <- 100
d <- 3

X <- matrix(runif(n * d, 0, 0.5), nrow = n, ncol = d)
X <- t(apply(X, 1, function(x) x / sqrt(sum(x^2)) * 0.7))

P_true <- X %*% t(X)
P_true <- pmin(pmax(P_true, 0), 1)
diag(P_true) <- 0

A <- matrix(0, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    if (runif(1) < P_true[i, j]) {
      A[i, j] <- 1
      A[j, i] <- 1
    }
  }
}

eigen_result <- eigen(A, symmetric = TRUE)
eigenvalues <- eigen_result$values[1:d]
eigenvectors <- eigen_result$vectors[, 1:d]

Lambda_sqrt <- diag(sqrt(pmax(eigenvalues, 0)))
X_hat <- eigenvectors %*% Lambda_sqrt

P_hat <- X_hat %*% t(X_hat)
P_hat <- pmin(pmax(P_hat, 0), 1)
diag(P_hat) <- 0

upper_tri_idx <- upper.tri(P_true)
P_true_vec <- P_true[upper_tri_idx]
P_hat_vec <- P_hat[upper_tri_idx]

MSE <- mean((P_true_vec - P_hat_vec)^2)
RMSE <- sqrt(MSE)
MAE <- mean(abs(P_true_vec - P_hat_vec))
correlation <- cor(P_true_vec, P_hat_vec)

cat(sprintf("MSE: %.6f\n", MSE))
cat(sprintf("RMSE: %.6f\n", RMSE))
cat(sprintf("MAE: %.6f\n", MAE))
cat(sprintf("Correlation: %.4f\n", correlation))

results <- list(
  X_true = X,
  X_hat = X_hat,
  P_true = P_true,
  P_hat = P_hat,
  A = A,
  MSE = MSE,
  RMSE = RMSE,
  MAE = MAE,
  correlation = correlation
)
