get_A <- function(beta_indices, length_beta) {
  # Initialize constraint matrix components
  non_neg_constraints <- list()
  monotonicity_constraints <- list()

  # 0. Reference category constraint (intercept set to 0)
  ref_constraint <- rep(0, length_beta)
  ref_constraint[1] <- 1  # Assuming the intercept or reference is at the first index
  non_neg_constraints <- append(non_neg_constraints, list(ref_constraint))

  # Set up grouped monotonicity constraints
  for (group in beta_indices) {
    n_levels <- length(group)

    # 1. Non-negativity constraint for the first level of this group
    #constraint <- rep(0, length_beta)
    #constraint[group[1]] <- 1
    #non_neg_constraints <- append(non_neg_constraints, list(constraint))

    # 2. Monotonicity constraints within the group
    for (k in 2:n_levels) {
      constraint <- rep(0, length_beta)
      constraint[group[k - 1]] <- -1
      constraint[group[k]] <- 1
      monotonicity_constraints <- append(monotonicity_constraints, list(constraint))
    }
  }

  # Combine all constraints
  Amat <- do.call(rbind, c(non_neg_constraints, monotonicity_constraints))

  return(Amat)
}

get_A <- function(nlevels_vec){
  monotonic <- list()
  prev_levels <- 0

  for(nlevel in nlevels_vec){
    constraint <- matrix(0, nrow = nlevel - 1, ncol = sum(nlevels_vec))
    for(i in 1:(nlevel - 1)){
      if(i == 1 & prev_levels > 0){
        constraint[i, 1] <- - 1
        constraint[i, i + prev_levels] <- 1
      } else{
        constraint[i, i + prev_levels - 1] <- -1
        constraint[i, i + prev_levels] <- 1
      }
    }
    prev_levels <- prev_levels + nlevel
    monotonic <- append(monotonic, list(constraint))
  }
  monotonic <- do.call(rbind, monotonic)
  return(monotonic)

}

get_A <- function(nlevels_vec) {
  # List to hold the monotonic constraint matrices
  monotonic <- list()

  # The total number of coefficients (including the intercept)
  total_levels <- sum(nlevels_vec)  # Total number of coefficients including the intercept

  prev_levels <- 0  # To track the previous levels for each variable

  # Loop over each categorical variable's levels
  for (nlevel in nlevels_vec) {
    # Create a constraint matrix for the current categorical variable (excluding intercept)
    constraint <- matrix(0, nrow = nlevel - 2, ncol = total_levels - (length(nlevels_vec) - 1))  # Total columns minus 1 (intercept)

    for (i in 1:(nlevel - 2)) {
      # For the first constraint, compare the first level to the second
      if (i == 1) {
        constraint[i, prev_levels + 2] <- -1  # No constraint on the intercept
        constraint[i, prev_levels + 3] <- 1  # First level (no constraint on intercept)
      } else {
        # For subsequent constraints, compare adjacent levels
        constraint[i, prev_levels + i + 1] <- -1  # Previous level
        constraint[i, prev_levels + i + 2] <- 1  # Current level
      }
    }

    # Update the previous level count for the next variable
    prev_levels <- prev_levels + (nlevel - 1)  # Exclude intercept from this count
    monotonic <- append(monotonic, list(constraint))  # Add the constraint matrix for the variable
  }

  # Combine the constraint matrices for all variables
  monotonic <- do.call(rbind, monotonic)

  # Ensure we have exactly 10 columns (1 for intercept, the rest for the variables)
  # Intercept is in the first column, and the non-intercept columns for the levels of the categorical variables
  return(monotonic)
}







#penalty function
penalty_vec_simp <- function(nlevels_vec){
  penalty <- c()
  for(nlevel in nlevels_vec){
    penalty <- c(penalty, rep(0, nlevel - 2), 1)
  }
  return(c(-1,penalty))
}

# Main function to set up and fit the penalized GLM model

strat_poisson <- function(df, y, vars, lambda=0, tol = 1e-6, max_iter = 100) {
  # Convert specified columns to factors
  nlevels <- c()
  for(var in vars) {
    df[[var]] <- as.factor(df[[var]])
    nlevels <- c(nlevels, length(levels(df[[var]])))
  }
  #df[[var1]] <- as.factor(df[[var1]])
  #df[[var2]] <- as.factor(df[[var2]])

  # Find the number of levels in var1 and var2
  #nlevels1 <- length(levels(df[[var1]]))
  #nlevels2 <- length(levels(df[[var2]]))

  # Initialize the design matrix, excluding the intercept
  #formula_str <- paste("~", var1, "+", var2, '-1')

  #Create design matrix for all vars - no intercept
  formula_str <- paste("~", paste(vars, collapse = " + "))

  X <- model.matrix(as.formula(formula_str), data = df)

  # Create beta_indices vector
  #var1_levels <- seq(nlevels1)  # Exclude reference category
  #var2_levels <- seq(nlevels1 + 1, nlevels1 + nlevels2)  # Adjusted for second factor

  #beta_indices <- list(var1_levels, var2_levels)

  # Total number of beta coefficients to estimate
  length_beta <- ncol(X)  # Use total columns in X for testing

  # Create the constraint matrix
  Amat <- get_A(nlevels)
  bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0

  y <- df[[y]]
  log_Y <- log(y + 0.1)

  # Initial coefficient estimates
  #beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y
  beta <- rep(1, ncol(X))
  epsilon <- 99
  iter <- 0

  #IRWLS
  while (epsilon > tol & iter <= max_iter) {
    eta <- X %*% beta
    mu <- exp(eta)
    nu <- exp(eta)

    A <- Matrix::.sparseDiagonal(as.vector(nu), n = nrow(nu))
    z <- eta + solve(A) %*% (y - mu)


    # Compute D and d for the quadratic program
    D <- 2 * t(X) %*% A %*% X  # Quadratic term
    d <- 2 * ((t(X) %*% A %*% z) - lambda *  penalty_vec_simp(nlevels)) # Linear term with penalty

    # Solve the quadratic program with constraints
    solution <- solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec)

    beta_new_star <- solution$solution  # Extract solution

    # Check convergence
    epsilon <- sqrt(t(beta_new_star - beta) %*% (beta_new_star - beta))
    beta <- beta_new_star  # Update beta

    #iterate
    iter <- iter + 1
    print(iter)
  }

  return(beta)
}

#test
res_simp <- strat_poisson(data, 'hospital_days', c('age_group', 'pulm_impair', 'num_comorb'), lambda = 0, tol = 1e-8, max_iter = 100)
length(res_simp)

## Get A function for model with interaction terms

# x <- model.matrix(~age_pulm, data)
#
# p <- length(levels(data$age_group))
# q <- length(levels(data$num_comorb))
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

#Constraint matrix function
get_Ax <- function(p, q) {
  nb <- 2*p*q - p - q #number of boundaries
  Ax <- matrix(0, nb, p*q)
  if(q > 1)
    for(i in 1:(p*(q - 1))) {
      Ax[i, i] <- -1
      Ax[i, i + p] <- 1
    }
  if(p > 1)
    for(i in 1:(q*(p - 1))) {
      j <- (i - 1)%/%(p-1) + i
      Ax[i + p*(q - 1), j:(j + 1)] <- c(-1, 1)
    }
  #Ax <- rbind(c(1, rep(0, p*q-1)), Ax)
  return(Ax)
}

#Penalty vector function
penalty_vec <- function(p, q){
  # Initialize a matrix with zeros
  mat <- matrix(0, nrow = p, ncol = q)

  # Assign values to the boundaries
  # Top boundary (first row, excluding the upper right corner)
  mat[1, ] <- -1

  # Left boundary (first column, excluding the top left corner)
  mat[-1, 1] <- -1

  # Bottom boundary (last row, excluding the bottom left corner)
  mat[p, ] <- 1

  # Right boundary (last column, excluding the upper right corner)
  mat[, q] <- 1

  #Fix corners
  mat[p, 1] <- 0
  mat[1, q] <- 0
  mat[p, q] <- 2

  # Convert the matrix to a vector
  boundary_vector <- as.vector(mat)
  return(boundary_vector)
}


#Main function
strat_poissonX <- function(df, y, var1, var2, lambda=0, tol = 1e-6, max_iter = 100, warm.start) {
  # Convert specified columns to factors
  df[[var1]] <- as.factor(df[[var1]])
  df[[var2]] <- as.factor(df[[var2]])

  # Find the number of levels in var1 and var2
  nlevels1 <- length(levels(df[[var1]]))
  nlevels2 <- length(levels(df[[var2]]))

  #print(paste0('p = ', nlevels1, ', q = ', nlevels2))

  #Create interaction terms
  df$interaction <- interaction(df[[var1]], df[[var2]])

  # Initialize the design matrix
  #X <- model.matrix(~interaction, data = df)
  X <- model.matrix(~interaction - 1, data = df) #cell means coding

  # Create the constraint matrix
  Amat <- get_Ax(p = nlevels1, q = nlevels2)
  bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0

  y <- df[[y]]

  #Find penalty vector
  penalty <- penalty_vec(p = nlevels1, q = nlevels2)

  # Initial coefficient estimates
  #beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y

  if (missing(warm.start)) {
    beta <- rep(1, ncol(X))
  } else {
    beta <- warm.start
  }

  epsilon <- 99
  iter <- 0

  #IRWLS
  while (epsilon > tol & iter <= max_iter) {
    eta <- X %*% beta
    mu <- exp(eta)
    nu <- exp(eta)

    A <- Matrix::.sparseDiagonal(as.vector(nu), n = nrow(nu))
    z <- eta + solve(A) %*% (y - mu)

    # Compute D and d for the quadratic program
    D <- 2 * t(X) %*% A %*% X  # Quadratic term
    #d <- 2 * t(X) %*% A %*% z + lambda * (beta[tail(beta_indices[[1]], 1)] + beta[tail(beta_indices[[2]], 1)])  # Linear term with penalty
    # d <- 2 * (t(X) %*% A %*% z - lambda * c(-1,-1, -1, -1, -1, -1, 0,
    #                                        -1, 0, 0, 0, 0, 0, 1,
    #                                        -1, 0, 0, 0, 0, 0, 1,
    #                                        0, 1, 1, 1, 1, 1, 2))
    d <- 2 * (t(X) %*% A %*% z - (lambda * penalty))

    # Solve the quadratic program with constraints
    solution <- solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec, meq = 0)

    beta_new_star <- solution$solution  # Extract solution

    # Check convergence
    epsilon <- sqrt(t(beta_new_star - beta) %*% (beta_new_star - beta))
    beta <- beta_new_star  # Update beta

    #iterate
    iter <- iter + 1
    #print(iter)
  }

  return(beta)
}


##Test out the functions
library(tidyverse)
hospital <- read_csv('covid19_hospitalizations.csv')
hospital %>%
  rename('hospital_days' = 'Hospitalization time in days',
         'age_group' = 'Age group',
         'num_comorb' = 'Total number of comorbidities',
         'risk_class' = 'Risk Classification Protocol',
         'pulm_impair' = 'Pulmonary impairment') %>%
  mutate(age_group = as.factor(age_group)) %>%
  mutate(num_comorb = as.factor(num_comorb)) %>%
  mutate(pulm_impair = as.factor(pulm_impair)) %>%
  mutate(age_pulm = interaction(age_group, pulm_impair)) %>%
  select(hospital_days, age_group, pulm_impair, age_pulm, num_comorb) %>%
  drop_na() -> data

levels(data$age_group) <- c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+')
levels(data$pulm_impair) <- c('None', 'Mild', 'Moderate', 'Severe')

resX <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = 0, tol = 1e-8, max_iter = 100)

matrix(resX, nrow = 7, ncol=4)

library(pheatmap)
pheatmap(matrix(exp(resX), nrow = 7, ncol=4), cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'),
         labels_col = c('None', 'Mild', 'Moderate', 'Severe'),
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status')



#Try CART
library(rpart)
library(rpart.plot)

# Fit a regression tree
cart_dat <- data
cart_dat$age_group <- as.ordered(cart_dat$age_group)
cart_dat$pulm_impair <- as.ordered(cart_dat$pulm_impair)
colnames(cart_dat) <- c('hospital_days', 'Age.Group', 'Pulmonary.Impairment', 'age_pulm')
fit <- rpart(hospital_days ~ Age.Group + Pulmonary.Impairment, data = cart_dat, method = "poisson")
summary(fit)
rpart.plot(fit)

#get predicted values for every combination of age group and pulmonary impairment
# Create a grid of all combinations of age group and pulmonary impairment
age_group_levels <- levels(data$age_group)
pulm_impair_levels <- levels(data$pulm_impair)
grid <- matrix(NA, nrow = length(age_group_levels), ncol = length(pulm_impair_levels))

for(age in age_group_levels) {
  for(pulm in pulm_impair_levels) {
    grid[age_group_levels == age, pulm_impair_levels == pulm] <- predict(fit, newdata = data.frame(age_group = age, pulm_impair = pulm))
  }
}

#Create heatmap
pheatmap(grid, cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = age_group_levels,
         labels_col = pulm_impair_levels,
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status')


#Create function to find BIC for a given model
bic_poisson <- function(df, y, var1, var2, beta) {
  df[[var1]] <- as.factor(df[[var1]])
  df[[var2]] <- as.factor(df[[var2]])
  y <- df[[y]]

  #Create interaction terms
  df$interaction <- interaction(df[[var1]], df[[var2]])

  # Initialize the design matrix
  #X <- model.matrix(~interaction, data = df)
  X <- model.matrix(~interaction - 1, data = df) #cell means coding
  n <- nrow(df)
  eta <- X %*% beta
  mu <- exp(eta)
  loglik <- sum(y * eta - mu - lfactorial(y))
  bic <- -2 * loglik + log(n) * length(unique(round(beta,8)))
  return(bic)
}


#Find optimal lambda over grid using BIC
bic_mse_lambda <- function(df, y, var1, var2, lambda_grid){
  bic_values <- c()
  n_unique_beta <- c()
  MSE <- c()

  df[[var1]] <- as.factor(df[[var1]])
  df[[var2]] <- as.factor(df[[var2]])
  y1 <- df[[y]]

  #Create interaction terms
  df$interaction <- interaction(df[[var1]], df[[var2]])

  # Initialize the design matrix with no intercept - cell means coding
  X <- model.matrix(~interaction - 1, data = df)

  for(i in 1:length(lambda_grid)){
    if(i == 1){
      beta <- strat_poissonX(df, y, var1, var2, lambda = lambda_grid[i], tol = 1e-8, max_iter = 100)
    }else{
      beta <- strat_poissonX(df, y, var1, var2, lambda = lambda_grid[i], tol = 1e-8, max_iter = 100, warm.start = beta)
    }

    #beta <- round(beta, 5)
    bic_values <- c(bic_values, bic_poisson(df, y, var1, var2, beta))
    n_unique_beta <- c(n_unique_beta, length(unique(beta)))
    MSE <- c(MSE, mean((exp(X %*% beta) - y1)^2))
  }
  result <- list(lambda.min = lambda_grid[which.min(bic_values)], bic_values = bic_values, MSE = MSE)

  return(result)
}

rez <- bic_mse_lambda(df=data, y = 'hospital_days', var1 = 'age_group', var2 = 'pulm_impair', lambda_grid=seq(from=0, to=10, by=0.01))
plot(x = seq(from=0, to=10, by=0.01), y = rez$bic_values, type='l')
lambda1 <- rez$lambda.min
lambda1

plot(x = seq(from=0, to=10, by=0.01), y = rez$MSE, type='l')

result_optim <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = lambda1, tol = 1e-8, max_iter = 100)
pheatmap(matrix(exp(result_optim), nrow = 7, ncol=4), cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'),
         labels_col = c('None', 'Mild', 'Moderate', 'Severe'),
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status',
         number_format = "%.2f")
other_result <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = 4.5, tol = 1e-8, max_iter = 100)
pheatmap(matrix(exp(other_result), nrow = 7, ncol=4), cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'),
         labels_col = c('None', 'Mild', 'Moderate', 'Severe'),
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status',
         number_format = "%.6f")

#Compare MSE between strat_poissonX vs CART
true <- data$hospital_days
#CART
pred_CART <- predict(fit)
cart_mse <- mean((pred_CART - true)^2)
cart_mse

#stratPoissonX
df <- data
y <- 'hospital_days'
var1 <- 'age_group'
var2 <- 'pulm_impair'

# Convert specified columns to factors
df[[var1]] <- as.factor(df[[var1]])
df[[var2]] <- as.factor(df[[var2]])

#Create interaction terms
df$interaction <- interaction(df[[var1]], df[[var2]])

# Initialize the design matrix with no intercept - cell means coding
X <- model.matrix(~interaction - 1, data = df)

#predictions
strat_pois_pred <- as.vector(exp(X %*% result_optim))

#mse
strat_pois_mse <- mean((strat_pois_pred - true)^2)
strat_pois_mse

#Create train/test data
set.seed(2024)
samp_size <- round(0.667 * nrow(data))
train_idx <- sample(1:nrow(data), size = samp_size, replace = F)
train <- data[train_idx, ]
test <- data[-train_idx, ]

#Train models

#CART
cart_fit <- rpart(hospital_days ~ age_group + pulm_impair, data = train, method = "poisson")

#stratPoissonX
pois_fit <- strat_poissonX(train, 'hospital_days', 'age_group', 'pulm_impair', lambda = lambda1, tol = 1e-8, max_iter = 100)

#find predicted vals
cart_pred <- predict(cart_fit, newdata = test)
#stratPoissonX
df <- test
y <- 'hospital_days'
var1 <- 'age_group'
var2 <- 'pulm_impair'

# Convert specified columns to factors
df[[var1]] <- as.factor(df[[var1]])
df[[var2]] <- as.factor(df[[var2]])

#Create interaction terms
df$interaction <- interaction(df[[var1]], df[[var2]])

# Initialize the design matrix with no intercept - cell means coding
X <- model.matrix(~interaction - 1, data = df)

#predictions
pois_pred <- as.vector(exp(X %*% pois_fit))

#MSE
true <- test$hospital_days
mean((cart_pred - true)^2)
mean((pois_pred - true)^2)

#evaluate run time
system.time(strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = lambda1, tol = 1e-8, max_iter = 100))
system.time(rpart(hospital_days ~ Age.Group + Pulmonary.Impairment, data = cart_dat, method = "poisson"))

#Effect of varying lambda
resl1 <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = 20, tol = 1e-8, max_iter = 100)

pheatmap(matrix(exp(resl1), nrow = 7, ncol=4), cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'),
         labels_col = c('None', 'Mild', 'Moderate', 'Severe'),
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status: Lambda = 20')

resl2 <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = 50, tol = 1e-8, max_iter = 100)

pheatmap(matrix(exp(resl2), nrow = 7, ncol=4), cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'),
         labels_col = c('None', 'Mild', 'Moderate', 'Severe'),
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status: Lambda = 50')

resl3 <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = 100, tol = 1e-8, max_iter = 100)

pheatmap(matrix(exp(resl3), nrow = 7, ncol=4), cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'),
         labels_col = c('None', 'Mild', 'Moderate', 'Severe'),
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status: Lambda = 100')


#ADMM attempt
library(quadprog)

# ADMM Algorithm
admm_poisson <- function(X, y, A, b, lambda, rho = 1, tol = 1e-6, max_iter = 200) {
  # Initialize variables
  p <- ncol(X)
  beta <- rep(1, p)
  z <- rep(1, p)
  u <- rep(1, p)

  for (k in 1:max_iter) {
    # Step 1: Update beta
    # Step 1: Update beta
    eta <- X %*% beta
    mu <- pmax(exp(eta), 1e-6)  # Bound mu to avoid instability
    W <- diag(as.vector(mu))    # Diagonal weight matrix
    Hessian <- t(X) %*% W %*% X + (lambda + rho) * diag(p)
    gradient <- -t(X) %*% (y - mu) + rho * (z - u)
    beta <- solve(Hessian, gradient)


    # Step 2: Update z
    z_update <- beta + u
    qp_solution <- solve.QP(Dmat = diag(p), dvec = z_update, Amat = t(A), bvec = b, meq = 0)
    z <- qp_solution$solution

    # Step 3: Update u
    u <- u + (beta - z)

    # Check convergence
    primal_residual <- sqrt(sum((beta - z)^2))
    dual_residual <- rho * sqrt(sum((z - z_update)^2))

    if (primal_residual < tol && dual_residual < tol) {
      cat("Converged at iteration:", k, "\n")
      break
    }
  }

  return(list(beta = beta, z = z, u = u))
}

# Example Usage
# X: Design matrix
# y: Observed counts
# A, b: Constraint matrix and vector
# lambda: Regularization parameter
# rho: ADMM penalty parameter

admm_poisson_debug <- function(X, y, A, b, lambda, rho = 0.1, tol = 1e-6, max_iter = 1000) {
  # Scale X
  X <- scale(X)

  # Initialize variables
  p <- ncol(X)
  beta <- solve(t(X) %*% X + lambda * diag(p)) %*% t(X) %*% log(y + 1)  # Ridge initialization
  z <- rep(0, p)
  u <- rep(0, p)

  for (k in 1:max_iter) {
    # Step 1: Update beta
    eta <- X %*% beta
    mu <- pmax(pmin(exp(eta), 1e6), 1e-6)  # Bound mu
    W <- diag(as.vector(mu))              # Diagonal weight matrix
    Hessian <- t(X) %*% W %*% X + (lambda + rho) * diag(p)
    gradient <- t(X) %*% (y - mu) + rho * (z - u)
    beta <- solve(Hessian, gradient)

    # Step 2: Update z
    z_update <- beta + u
    qp_solution <- solve.QP(Dmat = diag(p), dvec = z_update, Amat = t(A), bvec = b, meq = 0)
    z <- qp_solution$solution

    # Step 3: Update u
    u <- u + (beta - z)

    # Compute objective function and residuals for debugging
    obj <- sum(y * log(mu) - mu) - lambda * sum(beta^2)
    primal_residual <- sqrt(sum((beta - z)^2))
    dual_residual <- rho * sqrt(sum((z - z_update)^2))

    # Debugging outputs
    cat(sprintf("Iter: %d, Objective: %.6f, Primal Residual: %.6f, Dual Residual: %.6f\n",
                k, obj, primal_residual, dual_residual))
    cat("Beta Mean: ", mean(beta), " | Beta Range: ", range(beta), "\n")
    cat(sprintf("Iter: %d | Beta Norm: %.4f | Z Norm: %.4f | U Norm: %.4f\n",
                k, sqrt(sum(beta^2)), sqrt(sum(z^2)), sqrt(sum(u^2))))

    # Check convergence
    if (primal_residual < tol && dual_residual < tol) {
      cat("Converged at iteration:", k, "\n")
      break
    }
  }

  return(list(beta = beta, z = z, u = u))
}

admm_poisson_penalty <- function(X, y, A, b, lambda, rho, penalty=0.1, tol = 1e-6, max_iter = 100, verbose = TRUE) {
  p <- ncol(X)
  n <- nrow(X)

  # Initialize variables
  beta <- rep(0, p)  # Initial guess for beta
  z <- rep(0, p)     # Initial guess for z
  u <- rep(0, p)     # Dual variable

  iter <- 0
  epsilon_primal <- Inf
  epsilon_dual <- Inf

  while (iter < max_iter && (epsilon_primal > tol || epsilon_dual > tol)) {
    # Update beta
    # Update beta
    mu <- exp(X %*% beta)  # Predicted mean
    mu <- pmax(pmin(mu, 1e6), 1e-6)  # Clip mu to avoid extreme values
    W <- Matrix::.sparseDiagonal(x = as.vector(mu), n = nrow(X))  # Diagonal weight matrix

    # Add small ridge penalty for numerical stability
    epsilon <- 1e-6
    Hessian <- t(X) %*% W %*% X + rho * diag(p) + epsilon * diag(p)
    gradient <- t(X) %*% (y - mu) + rho * (z - u) - lambda * penalty

    # Solve the linear system
    beta <- solve(Hessian, gradient)


    # Update z (projection onto feasible set Az >= b)
    z_old <- z
    z <- pmax(beta + u, 0)  # Ensure non-negativity (or apply other projection)

    # Update dual variable u
    u <- u + beta - z

    # Calculate residuals
    epsilon_primal <- sqrt(sum((beta - z)^2))
    epsilon_dual <- rho * sqrt(sum((z - z_old)^2))

    # Print progress if verbose
    if (verbose) {
      objective <- sum(y * (X %*% beta) - exp(X %*% beta)) - lambda * sum(penalty * beta)
      cat(sprintf(
        "Iter: %d, Objective: %.6f, Primal Residual: %.6f, Dual Residual: %.6f\n",
        iter + 1, objective, epsilon_primal, epsilon_dual
      ))
      cat(sprintf("Beta Mean:  %.6f  | Beta Range:  %.6f %.6f \n", mean(beta), min(beta), max(beta)))
    }

    iter <- iter + 1
  }

  return(list(beta = beta, z = z, u = u, iter = iter))
}

admm_poisson <- function(X, y, A, b, lambda, rho, tol = 1e-4, max_iter = 100) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)  # Initialize beta
  z <- rep(0, p)     # Initialize z
  u <- rep(0, p)     # Initialize dual variable

  iter <- 1
  while (iter <= max_iter) {
    # Step 1: Update beta
    mu <- exp(X %*% beta)
    W <- Matrix::.sparseDiagonal(x = pmax(mu, 1e-6), n = n)
    eta <- X %*% beta + (y - mu) / pmax(mu, 1e-6)

    Hessian <- t(X) %*% W %*% X + rho * diag(p) + 1e-6 * diag(p)
    gradient <- -t(X) %*% (y - mu) + rho * (z - u) - lambda
    beta <- solve(Hessian, gradient)

    # Step 2: Update z (Projection onto constraints)
    z <- pmax(A %*% beta + u, b)

    # Step 3: Update u
    u <- u + beta - z

    # Check convergence
    primal_residual <- sqrt(sum((beta - z)^2))
    dual_residual <- sqrt(sum((z - z_old)^2))
    if (primal_residual < tol && dual_residual < tol) break

    iter <- iter + 1
  }
  return(list(beta = beta, iter = iter, z = z, u = u))
}


# Example dimensions
set.seed(42)
X <- model.matrix(~interaction - 1, data = df)
y <- df$hospital_days
nlevels1 <- length(levels(df$age_group))
nlevels2 <- length(levels(df$pulm_impair))
A <- get_Ax(p = nlevels1, q = nlevels2)
b <- rep(0, nrow(A))
lambda <- 0.1
rho <- 1

result <- admm_poisson(X, y, A, b, lambda = 0.1, rho=1/2)
print(result$beta)

exp(matrix(result$beta, nrow = 7, ncol=4))




