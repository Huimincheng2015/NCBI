library(ggplot2)

set.seed(123)

n <- 100
d <- 3
n_reps <- 30
hyper_values <- 1:10

results_df <- data.frame()

for (hyper in hyper_values) {
  mse_values <- numeric(n_reps)

  for (rep in 1:n_reps) {
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
    eigenvalues <- eigen_result$values[1:hyper]
    eigenvectors <- eigen_result$vectors[, 1:hyper]

    if (hyper == 1) {
      Lambda_sqrt <- sqrt(max(eigenvalues, 0))
      X_hat <- eigenvectors * Lambda_sqrt
    } else {
      Lambda_sqrt <- diag(sqrt(pmax(eigenvalues, 0)))
      X_hat <- eigenvectors %*% Lambda_sqrt
    }

    P_hat <- X_hat %*% t(X_hat)
    P_hat <- pmin(pmax(P_hat, 0), 1)
    diag(P_hat) <- 0

    upper_tri_idx <- upper.tri(P_true)
    P_true_vec <- P_true[upper_tri_idx]
    P_hat_vec <- P_hat[upper_tri_idx]

    mse_values[rep] <- mean((P_true_vec - P_hat_vec)^2)
  }

  avg_mse <- mean(mse_values)
  sd_mse <- sd(mse_values)

  results_df <- rbind(results_df, data.frame(
    hyper = hyper,
    avg_mse = avg_mse,
    sd_mse = sd_mse
  ))

  cat(sprintf("hyper = %2d: MSE = %.6f (SD = %.6f)\n", hyper, avg_mse, sd_mse))
}

p <- ggplot(results_df, aes(x = hyper, y = avg_mse)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  geom_vline(xintercept = d, color = "red", linetype = "dashed", size = 1) +
  geom_errorbar(aes(ymin = avg_mse - sd_mse, ymax = avg_mse + sd_mse),
                width = 0.2, alpha = 0.5) +
  labs(title = "MSE vs Hyperparameter (d = 3)",
       x = "Hyperparameter (dimension)",
       y = "Average MSE") +
  annotate("text", x = d, y = max(results_df$avg_mse),
           label = sprintf("True d = %d", d), color = "red", vjust = -0.5) +
  theme_minimal()

ggsave("mse_vs_hyper.png", p, width = 8, height = 6)
cat("\nPlot saved as mse_vs_hyper.png\n")

print(results_df)
