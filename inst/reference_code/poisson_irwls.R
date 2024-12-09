# Define the function to calculate the log-likelihood
# log_likelihood <- function(beta, X, y) {
#   lambda <- exp(X %*% beta)
#   ll <- sum(y * log(lambda) - lambda - lgamma(y + 1))
#   return(ll)  # Return the log-likelihood
# }
#
# # Define the BIC calculation function
# calculate_bic <- function(log_likelihood, beta, n) {
#   # Degrees of freedom k is the number of non-zero parameters
#   k <- sum(beta != 0)
#   BIC <- -2 * log_likelihood + k * log(n)
#   return(BIC)
# }
#
# # IRWLS iterative process
# IRWLS_poisson <- function(X, y, tol = 1e-6, max_iter = 100) {
#   p <- ncol(X)
#   beta <- rep(0, p)  # Initialize beta
#
#   for (iter in 1:max_iter) {
#     # Update weights and working variables
#     weights <- update_weights(beta, X, y)
#     W <- weights$W
#     z <- weights$z
#
#     # Use weighted least squares to update beta
#     beta_new <- solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% z)
#
#     # Check for convergence
#     if (sum(abs(beta_new - beta)) < tol) {
#       break
#     }
#
#     beta <- beta_new  # Update beta
#   }
#
#   return(beta)
# }
#
# # IRWLS_poisson_pen <- function(X, y, lambda, tol = 1e-6, max_iter = 100) {
# #   n=nrow(df)
# #   log_Y=log(y+0.1)
# #   q=ncol(X)
# #   beta=solve(t(X) %*% X) %*% t(X) %*% log_Y
# #   epsilon=99
# #   iter=0
# #   while (epsilon > tol & iter <= max_iter){
# #     eta=X %*% beta
# #     mu=exp(eta)
# #     nu=exp(eta)
# #     V=diag(x=as.vector(nu))
# #     Z=eta+solve(V) %*% (y-mu)
# #     # Compute D and d for the quadratic program
# #     D <- t(X) %*% V %*% X                   # Quadratic term
# #     d <- t(X) %*% V %*% Z - lambda * rep(1, q)  # Linear term with penalty
# #
# #     # Solve the quadratic program
# #     beta_new_star <- solve.QP(Dmat = D, dvec = d, Amat = NULL, bvec = NULL)
# #
# #     #beta_new=solve(t(X) %*% V %*% X) %*% t(X) %*% V %*% Z
# #     epsilon = sqrt(t(beta_new_star-beta)%*%(beta_new_star-beta))
# #     beta=beta_new_star
# #     iter=iter+1
# #     beta_t=t(beta)
# #     print(paste("ite:",iter,"beta_0:", beta[1], "beta_1:", beta[2], "beta_2:", beta[3], "beta_3:",
# #                 beta[4],"epsilon:", epsilon))
# #   }
# #
# #   return(beta)
# # }
#
# IRWLS_poisson_pen <- function(X, y, lambda, tol = 1e-6, max_iter = 100) {
#   n <- nrow(X)
#   q <- ncol(X)  # Total number of covariate-level combinations
#   log_Y <- log(y + 0.1)
#
#   # Initial coefficient estimates
#   beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y
#   epsilon <- 99
#   iter <- 0
#
#   # Construct Amat and bvec for constraints
#   # 1. Non-negativity constraints: beta >= 0
#   non_neg <- diag(q)
#
#   # 2. Monotonicity constraints: beta_j_k - beta_j_(k-1) >= 0
#   monotonicity <- matrix(0, nrow = q - 1, ncol = q)
#   for (i in 1:(q - 1)) {
#     monotonicity[i, i] <- -1
#     monotonicity[i, i + 1] <- 1
#   }
#
#   # Combine constraints into Amat
#   Amat <- rbind(non_neg, monotonicity)
#   bvec <- rep(0, nrow(Amat))  # All constraints are >= 0
#
#   while (epsilon > tol & iter <= max_iter) {
#     eta <- X %*% beta
#     mu <- exp(eta)
#     nu <- exp(eta)
#
#     V <- diag(as.vector(nu))
#     Z <- eta + solve(V) %*% (y - mu)
#
#     # Compute D and d for the quadratic program
#     D <- t(X) %*% V %*% X  # Quadratic term
#     d <- t(X) %*% V %*% Z - lambda * rep(1, q)  # Linear term with penalty
#
#     # Solve the quadratic program with constraints
#     solution <- solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec, meq = 0)
#
#     beta_new_star <- solution$solution  # Extract solution
#
#     # Check convergence
#     epsilon <- sqrt(t(beta_new_star - beta) %*% (beta_new_star - beta))
#     beta <- beta_new_star  # Update beta
#
#     # Print progress
#     print(paste(
#       "iter:", iter,
#       "beta_0:", beta[1], "beta_1:", beta[2],
#       "epsilon:", epsilon
#     ))
#
#     iter <- iter + 1
#   }
#
#   return(beta)
# }
#
#
# n=nrow(df)
# Y=df$Y
# log_Y=log(df$Y+0.1)
# X=as.matrix(cbind(1,df[,c("age","base","Z")]))
# q=ncol(X)
# beta=solve(t(X) %*% X) %*% t(X) %*% log_Y
# tol=0.00001
# epsilon=99
# ite_max=25
# ite=0
# while (epsilon > tol & ite <= ite_max){
#   eta=X %*% beta
#   mu=exp(eta)
#   nu=exp(eta)
#   V=diag(x=as.vector(nu))
#   Z=eta+solve(V) %*% (Y-mu)
#   # Compute D and d for the quadratic program
#   D <- t(X) %*% V %*% X                   # Quadratic term
#   d <- t(X) %*% V %*% Z - lambda * rep(1, q)  # Linear term with penalty
#
#   # Solve the quadratic program
#   beta_new_star <- solve.QP(Dmat = D, dvec = d, Amat = NULL, bvec = NULL)
#
#   beta_new=solve(t(X) %*% V %*% X) %*% t(X) %*% V %*% Z
#   epsilon = sqrt(t(beta_new_star-beta)%*%(beta_new_star-beta))
#   beta=beta_new_star
#   ite=ite+1
#   beta_t=t(beta)
#   print(paste("ite:",ite,"beta_0:", beta[1], "beta_1:", beta[2], "beta_2:", beta[3], "beta_3:",
#               beta[4],"epsilon:", epsilon))
# }
#
# # Main program: Regularization parameter selection using BIC
# poisson_with_bic <- function(X, y, lambda_seq, tol = 1e-6, max_iter = 100) {
#   n <- nrow(X)  # Number of samples
#   best_bic <- Inf  # Initialize best BIC
#   best_lambda <- NULL  # Best lambda
#   best_beta <- NULL  # Best beta
#
#   for (lambda in lambda_seq) {
#     # Fit the model using IRWLS
#     beta <- IRWLS_poisson_pen(X, y, lambda, tol ,max_iter)
#
#     # Calculate the log-likelihood
#     log_lik <- log_likelihood(beta, X, y)
#
#     # Calculate BIC
#     bic <- calculate_bic(log_lik, beta, n)
#
#     # If the current BIC is smaller, update the best result
#     if (bic < best_bic) {
#       best_bic <- bic
#       best_lambda <- lambda
#       best_beta <- beta
#     }
#   }
#
#   # Return the best result
#   return(list(best_lambda = best_lambda, best_beta = best_beta, best_bic = best_bic))
# }
#
# # Example data
# set.seed(42)
# n <- 100
# p <- 3
# X <- matrix(rnorm(n * p), nrow = n, ncol = p)
# beta_true <- c(-1, 1, 1.5)
# y <- rpois(n, lambda = round(exp(X %*% beta_true)))
# IRWLS_poisson_pen(X, y, lambda=0.2)
#
# X <- rnorm(100)
# beta_true <-
#
# # Define a sequence of lambda values to search over
# lambda_seq <- seq(0.01, 1, length.out = 10)
#
# # Use BIC to select the best lambda
# result <- poisson_with_bic(X, y, lambda_seq)
# print(result$best_lambda)
# print(result$best_beta)
# print(result$best_bic)
