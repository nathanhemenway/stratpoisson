#test code here

# library(NHANES)
# data("NHANES")
# is.ordered(NHANES$BMI_WHO)
# NHANES$BMI_WHO <- factor(NHANES$BMI_WHO, ordered = TRUE)
# NHANES$BMI_WHO
#
#
#
# #Chat GPT update:
#
# # Function to build the constraint matrix A:
#
# get_A <- function(beta_indices, length_beta) {
#   # Initialize constraint matrix components
#   non_neg_constraints <- list()
#   monotonicity_constraints <- list()
#
#   # 0. Reference category constraint (intercept set to 0)
#   ref_constraint <- rep(0, length_beta)
#   ref_constraint[1] <- 1  # Assuming the intercept or reference is at the first index
#   non_neg_constraints <- append(non_neg_constraints, list(ref_constraint))
#
#   # Set up grouped monotonicity constraints
#   for (group in beta_indices) {
#     n_levels <- length(group)
#
#     # 1. Non-negativity constraint for the first level of this group
#     constraint <- rep(0, length_beta)
#     constraint[group[1]] <- 1
#     non_neg_constraints <- append(non_neg_constraints, list(constraint))
#
#     # 2. Monotonicity constraints within the group
#     for (k in 2:n_levels) {
#       constraint <- rep(0, length_beta)
#       constraint[group[k - 1]] <- -1
#       constraint[group[k]] <- 1
#       monotonicity_constraints <- append(monotonicity_constraints, list(constraint))
#     }
#   }
#
#   # Combine all constraints
#   Amat <- do.call(rbind, c(non_neg_constraints, monotonicity_constraints))
#
#   return(Amat)
# }
#
# # Main function to set up and fit the penalized GLM model
#
# strat_poisson <- function(df, y, var1, var2, lambda=0, tol = 1e-6, max_iter = 100) {
#   # Convert specified columns to factors
#   df[[var1]] <- as.factor(df[[var1]])
#   df[[var2]] <- as.factor(df[[var2]])
#
#   # Find the number of levels in var1 and var2
#   nlevels1 <- length(levels(df[[var1]]))
#   nlevels2 <- length(levels(df[[var2]]))
#
#   # Initialize the design matrix, excluding the intercept
#   formula_str <- paste("~", var1, "+", var2)
#   X <- model.matrix(as.formula(formula_str), data = df)
#
#   # Create beta_indices vector: identifies which columns in X correspond to each non-reference level
#   var1_levels <- seq(2, nlevels1)  # Exclude reference category
#   var2_levels <- seq(nlevels1 + 1, nlevels1 + nlevels2 - 1)  # Adjusted for second factor
#
#   beta_indices <- list(var1_levels, var2_levels)
#
#   # Total number of beta coefficients to estimate
#   length_beta <- ncol(X)  # Use total columns in X for testing
#
#   # Create the constraint matrix
#   Amat <- get_A(beta_indices, length_beta)
#   bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0
#
#   y <- df[[y]]
#   log_Y <- log(y + 0.1)
#
#   # Initial coefficient estimates
#   #beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y
#   beta <- c(0, rep(1, ncol(X)-1))
#   epsilon <- 99
#   iter <- 0
#
#     #IRWLS
#     while (epsilon > tol & iter <= max_iter) {
#       eta <- X %*% beta
#       mu <- exp(eta)
#       nu <- exp(eta)
#
#       A <- .sparseDiagonal(as.vector(nu), n = nrow(nu))
#       z <- eta + solve(A) %*% (y - mu)
#
#       # Compute D and d for the quadratic program
#       D <- 2 * t(X) %*% A %*% X  # Quadratic term
#       d <- 2 * t(X) %*% A %*% z - lambda * (beta[tail(beta_indices[[1]], 1)] + beta[tail(beta_indices[[2]], 1)])  # Linear term with penalty
#
#       # Solve the quadratic program with constraints
#       solution <- solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec, meq = 1)
#
#       beta_new_star <- solution$solution  # Extract solution
#
#       # Check convergence
#       epsilon <- sqrt(t(beta_new_star - beta) %*% (beta_new_star - beta))
#       beta <- beta_new_star  # Update beta
#
#       #iterate
#       iter <- iter + 1
#       print(iter)
#     }
#
#     return(beta)
# }
#
# #Test out the function
# strat_poisson(df, 'SleepHrsNight', 'BMI_WHO', 'AgeDecade', lambda = 0.1, tol = 1e-6, max_iter = 10)
#
# library(NHANES)
# data('NHANES')
#
# library(tidyverse)
# NHANES %>%
#   select(SleepHrsNight, BMI_WHO, AgeDecade) %>%
#   drop_na() -> df
#
# lambda <- 0.1
# var1 <- 'BMI_WHO'
# var2 <- 'AgeDecade'
# df[[var1]] <- as.factor(df[[var1]])
# df[[var2]] <- as.factor(df[[var2]])
#
# # Find the number of levels in var1 and var2
# nlevels1 <- length(levels(df[[var1]]))
# nlevels2 <- length(levels(df[[var2]]))
#
# # Initialize the design matrix, excluding the intercept
# formula_str <- paste("~", var1, "+", var2)
# X <- model.matrix(as.formula(formula_str), data = df)
# y <- df[['SleepHrsNight']]
# log_Y <- log(y + 1)
#
# # Create beta_indices vector: identifies which columns in X correspond to each non-reference level
# var1_levels <- seq(2, nlevels1)  # Exclude reference category
# var2_levels <- seq(nlevels1 + 1, nlevels1 + nlevels2 - 1)  # Adjusted for second factor
#
# beta_indices <- list(var1_levels, var2_levels)
#
# # Total number of beta coefficients to estimate
# length_beta <- ncol(X)  # Use total columns in X for testing
#
# # Create the constraint matrix
# Amat <- get_A(beta_indices, length_beta)
# bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0
#
# beta <- c(0, rep(1, ncol(X)-1))
# eta <- X %*% beta
# mu <- exp(eta)
# nu <- exp(eta)
#
# A <- .sparseDiagonal(as.vector(nu), n = nrow(nu))
# z <- eta + solve(A) %*% (y - mu)
#
# # Compute D and d for the quadratic program
# D <- 2 * t(X) %*% A %*% X  # Quadratic term
# d <- 2 * t(X) %*% A %*% z - lambda * (beta[tail(beta_indices[[1]], 1)] + beta[tail(beta_indices[[2]], 1)])  # Linear term with penalty
#
# # Solve the quadratic program with constraints
# library(quadprog)
# solution <- solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec, meq = 1)
# solution$solution
#
#
# # Initial coefficient estimates
# beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y
# dim(X)
#
# print(c(dim(X), dim(log_Y)))
#
#
# ##Try running on hospitalization data
# hospital <- read_csv('covid19_hospitalizations.csv')
#
# hospital %>%
#   select('Hospitalization time in days', 'Age group', 'Total number of comorbidities') %>%
#   rename('hospital_days' = 'Hospitalization time in days',
#          'age_group' = 'Age group',
#          'num_comorb' = 'Total number of comorbidities') %>%
#   drop_na() -> hospital
#
# hospital_res <- strat_poisson(hospital, 'hospital_days', 'age_group', 'num_comorb', lambda=0.1, tol = 1e-6, max_iter = 100)
#
# #make table to look at results
# age_effects <- hospital_res[1:7]
# comorb_effects <- hospital_res[8:14]
#
# table_results <- data.frame(age_group = 1:7,
#                              age_effect = age_effects)
# library(kableExtra)
# table_results %>%
#   kable(digits = 2) %>%
#   kable_styling(full_width = F)
#
# #look at comorb results
# table_results <- data.frame(num_comorb = 0:7,
#                              comorb_effect = c(0,comorb_effects))
#
# table_results %>%
#   kable(digits = 2) %>%
#   kable_styling(full_width = F)
#
#
# #Look at interaction matrix function from GLIDARS
# p <- 3
# q <- 3
# nb <- 2*p*q - p - q
# Apo0 <- matrix(0, nb, p*q)
# if(q > 1)
#   for(i in 1:(p*(q - 1))) {
#     Apo0[i, i] <- -1
#     Apo0[i, i + p] <- 1
#   }
# if(p > 1)
#   for(i in 1:(q*(p - 1))) {
#     j <- (i - 1)%/%(p-1) + i
#     Apo0[i + p*(q - 1), j:(j + 1)] <- c(-1, 1)
#   }
#
#
# #Try to make correct penalty vector
# nrow <- 7
# ncol <- 4
#
# # Initialize a matrix with zeros
# mat <- matrix(0, nrow = nrow, ncol = ncol)
#
# # Assign values to the boundaries
# # Top boundary (first row, excluding the upper right corner)
# mat[1, ] <- -1
#
# # Left boundary (first column, excluding the top left corner)
# mat[-1, 1] <- -1
#
# # Bottom boundary (last row, excluding the bottom left corner)
# mat[nrow, ] <- 1
#
# # Right boundary (last column, excluding the upper right corner)
# mat[, ncol] <- 1
#
# #Fix corners
# mat[nrow, 1] <- 0
# mat[1, ncol] <- 0
#
# # Convert the matrix to a vector
# boundary_vector <- as.vector(mat)
#
#
# #example run through function
# df <- data
# y <- 'hospital_days'
# var1 <- 'age_group'
# var2 <- 'pulm_impair'
#
# # Convert specified columns to factors
# df[[var1]] <- as.factor(df[[var1]])
# df[[var2]] <- as.factor(df[[var2]])
#
# # Find the number of levels in var1 and var2
# nlevels1 <- length(levels(df[[var1]]))
# nlevels2 <- length(levels(df[[var2]]))
#
# print(paste0('p = ', nlevels1, ', q = ', nlevels2))
#
# #Create interaction terms
# df$interaction <- interaction(df[[var1]], df[[var2]])
#
# # Initialize the design matrix with no intercept - cell means coding
# X <- model.matrix(~interaction - 1, data = df)
# X <- model.matrix(~interaction, data = df)
#
