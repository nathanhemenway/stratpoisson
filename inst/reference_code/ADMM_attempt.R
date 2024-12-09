#ADMM

#ADMM attempt
# library(quadprog)
#
# # ADMM Algorithm
# admm_poisson <- function(X, y, A, b, lambda, rho = 1, tol = 1e-6, max_iter = 200) {
#   # Initialize variables
#   p <- ncol(X)
#   beta <- rep(1, p)
#   z <- rep(1, p)
#   u <- rep(1, p)
#
#   for (k in 1:max_iter) {
#     # Step 1: Update beta
#     # Step 1: Update beta
#     eta <- X %*% beta
#     mu <- pmax(exp(eta), 1e-6)  # Bound mu to avoid instability
#     W <- diag(as.vector(mu))    # Diagonal weight matrix
#     Hessian <- t(X) %*% W %*% X + (lambda + rho) * diag(p)
#     gradient <- -t(X) %*% (y - mu) + rho * (z - u)
#     beta <- solve(Hessian, gradient)
#
#
#     # Step 2: Update z
#     z_update <- beta + u
#     qp_solution <- solve.QP(Dmat = diag(p), dvec = z_update, Amat = t(A), bvec = b, meq = 0)
#     z <- qp_solution$solution
#
#     # Step 3: Update u
#     u <- u + (beta - z)
#
#     # Check convergence
#     primal_residual <- sqrt(sum((beta - z)^2))
#     dual_residual <- rho * sqrt(sum((z - z_update)^2))
#
#     if (primal_residual < tol && dual_residual < tol) {
#       cat("Converged at iteration:", k, "\n")
#       break
#     }
#   }
#
#   return(list(beta = beta, z = z, u = u))
# }
#
# # Example Usage
# # X: Design matrix
# # y: Observed counts
# # A, b: Constraint matrix and vector
# # lambda: Regularization parameter
# # rho: ADMM penalty parameter
#
# admm_poisson_debug <- function(X, y, A, b, lambda, rho = 0.1, tol = 1e-6, max_iter = 1000) {
#   # Scale X
#   X <- scale(X)
#
#   # Initialize variables
#   p <- ncol(X)
#   beta <- solve(t(X) %*% X + lambda * diag(p)) %*% t(X) %*% log(y + 1)  # Ridge initialization
#   z <- rep(0, p)
#   u <- rep(0, p)
#
#   for (k in 1:max_iter) {
#     # Step 1: Update beta
#     eta <- X %*% beta
#     mu <- pmax(pmin(exp(eta), 1e6), 1e-6)  # Bound mu
#     W <- diag(as.vector(mu))              # Diagonal weight matrix
#     Hessian <- t(X) %*% W %*% X + (lambda + rho) * diag(p)
#     gradient <- t(X) %*% (y - mu) + rho * (z - u)
#     beta <- solve(Hessian, gradient)
#
#     # Step 2: Update z
#     z_update <- beta + u
#     qp_solution <- solve.QP(Dmat = diag(p), dvec = z_update, Amat = t(A), bvec = b, meq = 0)
#     z <- qp_solution$solution
#
#     # Step 3: Update u
#     u <- u + (beta - z)
#
#     # Compute objective function and residuals for debugging
#     obj <- sum(y * log(mu) - mu) - lambda * sum(beta^2)
#     primal_residual <- sqrt(sum((beta - z)^2))
#     dual_residual <- rho * sqrt(sum((z - z_update)^2))
#
#     # Debugging outputs
#     cat(sprintf("Iter: %d, Objective: %.6f, Primal Residual: %.6f, Dual Residual: %.6f\n",
#                 k, obj, primal_residual, dual_residual))
#     cat("Beta Mean: ", mean(beta), " | Beta Range: ", range(beta), "\n")
#     cat(sprintf("Iter: %d | Beta Norm: %.4f | Z Norm: %.4f | U Norm: %.4f\n",
#                 k, sqrt(sum(beta^2)), sqrt(sum(z^2)), sqrt(sum(u^2))))
#
#     # Check convergence
#     if (primal_residual < tol && dual_residual < tol) {
#       cat("Converged at iteration:", k, "\n")
#       break
#     }
#   }
#
#   return(list(beta = beta, z = z, u = u))
# }
#
# admm_poisson_penalty <- function(X, y, A, b, lambda, rho, penalty=0.1, tol = 1e-6, max_iter = 100, verbose = TRUE) {
#   p <- ncol(X)
#   n <- nrow(X)
#
#   # Initialize variables
#   beta <- rep(0, p)  # Initial guess for beta
#   z <- rep(0, p)     # Initial guess for z
#   u <- rep(0, p)     # Dual variable
#
#   iter <- 0
#   epsilon_primal <- Inf
#   epsilon_dual <- Inf
#
#   while (iter < max_iter && (epsilon_primal > tol || epsilon_dual > tol)) {
#     # Update beta
#     # Update beta
#     mu <- exp(X %*% beta)  # Predicted mean
#     mu <- pmax(pmin(mu, 1e6), 1e-6)  # Clip mu to avoid extreme values
#     W <- Matrix::.sparseDiagonal(x = as.vector(mu), n = nrow(X))  # Diagonal weight matrix
#
#     # Add small ridge penalty for numerical stability
#     epsilon <- 1e-6
#     Hessian <- t(X) %*% W %*% X + rho * diag(p) + epsilon * diag(p)
#     gradient <- t(X) %*% (y - mu) + rho * (z - u) - lambda * penalty
#
#     # Solve the linear system
#     beta <- solve(Hessian, gradient)
#
#
#     # Update z (projection onto feasible set Az >= b)
#     z_old <- z
#     z <- pmax(beta + u, 0)  # Ensure non-negativity (or apply other projection)
#
#     # Update dual variable u
#     u <- u + beta - z
#
#     # Calculate residuals
#     epsilon_primal <- sqrt(sum((beta - z)^2))
#     epsilon_dual <- rho * sqrt(sum((z - z_old)^2))
#
#     # Print progress if verbose
#     if (verbose) {
#       objective <- sum(y * (X %*% beta) - exp(X %*% beta)) - lambda * sum(penalty * beta)
#       cat(sprintf(
#         "Iter: %d, Objective: %.6f, Primal Residual: %.6f, Dual Residual: %.6f\n",
#         iter + 1, objective, epsilon_primal, epsilon_dual
#       ))
#       cat(sprintf("Beta Mean:  %.6f  | Beta Range:  %.6f %.6f \n", mean(beta), min(beta), max(beta)))
#     }
#
#     iter <- iter + 1
#   }
#
#   return(list(beta = beta, z = z, u = u, iter = iter))
# }
#
# admm_poisson <- function(X, y, A, b, lambda, rho, tol = 1e-4, max_iter = 100) {
#   n <- nrow(X)
#   p <- ncol(X)
#   beta <- rep(0, p)  # Initialize beta
#   z <- rep(0, p)     # Initialize z
#   u <- rep(0, p)     # Initialize dual variable
#
#   iter <- 1
#   while (iter <= max_iter) {
#     # Step 1: Update beta
#     mu <- exp(X %*% beta)
#     W <- Matrix::.sparseDiagonal(x = pmax(mu, 1e-6), n = n)
#     eta <- X %*% beta + (y - mu) / pmax(mu, 1e-6)
#
#     Hessian <- t(X) %*% W %*% X + rho * diag(p) + 1e-6 * diag(p)
#     gradient <- -t(X) %*% (y - mu) + rho * (z - u) - lambda
#     beta <- solve(Hessian, gradient)
#
#     # Step 2: Update z (Projection onto constraints)
#     z <- pmax(A %*% beta + u, b)
#
#     # Step 3: Update u
#     u <- u + beta - z
#
#     # Check convergence
#     primal_residual <- sqrt(sum((beta - z)^2))
#     dual_residual <- sqrt(sum((z - z_old)^2))
#     if (primal_residual < tol && dual_residual < tol) break
#
#     iter <- iter + 1
#   }
#   return(list(beta = beta, iter = iter, z = z, u = u))
# }
#
#
# # Example dimensions
# set.seed(42)
# X <- model.matrix(~interaction - 1, data = df)
# y <- df$hospital_days
# nlevels1 <- length(levels(df$age_group))
# nlevels2 <- length(levels(df$pulm_impair))
# A <- get_Ax(p = nlevels1, q = nlevels2)
# b <- rep(0, nrow(A))
# lambda <- 0.1
# rho <- 1
#
# result <- admm_poisson(X, y, A, b, lambda = 0.1, rho=1/2)
# print(result$beta)
#
# exp(matrix(result$beta, nrow = 7, ncol=4))
#
#
#
#
