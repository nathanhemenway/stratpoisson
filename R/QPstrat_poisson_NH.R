
#' get_A
#'
#' @description Creates the matrix of constraints to be passed into quadratic programming problem for non-interaction model
#' @param nlevels_vec This is an integer vector containing the number of levels for each categorical variable
#' @return An n x p matrix of constraints, where n is the number of constraints, and p is the length of beta
#' @examples
#' get_A(c(7, 4))
#' @export
get_A <- function(nlevels_vec) {
  # List to hold the monotonic constraint matrices
  monotonic <- list()

  # The total number of coefficients (including the intercept)
  total_levels <- sum(nlevels_vec)  # Total number of coefficients including the intercept

  prev_levels <- 0  # To track the previous levels for each variable

  # Loop over each categorical variable's levels
  for (nlevel in nlevels_vec) {
    # Create a constraint matrix for the current categorical variable (excluding intercept)
    constraint <- matrix(0, nrow = nlevel - 1, ncol = total_levels - (length(nlevels_vec) - 1))  # Total columns minus 1 (intercept)

    for (i in 1:(nlevel - 1)) {
      # For the first constraint, compare the first level to the second
      if (i == 1) {
        constraint[i, prev_levels + 1] <- 0  # No constraint on the intercept
        constraint[i, prev_levels + 2] <- 1  # First level (no constraint on intercept)
      } else {
        # For subsequent constraints, compare adjacent levels
        constraint[i, prev_levels + i] <- -1  # Previous level
        constraint[i, prev_levels + i + 1] <- 1  # Current level
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



#' penalty_vec_simp
#'
#' @description Creates the penalty vector for the non-interaction model
#' @param nlevels_vec This is an integer vector containing the number of levels for each categorical variable
#' @return a vector of length beta that when multiplied by beta gives the LASSO penalty
#' @examples
#' penalty_vec_simp(c(7,4))
#' @export
#'
#penalty function for non-interaction
penalty_vec_simp <- function(nlevels_vec){
  penalty <- c()
  for(nlevel in nlevels_vec){
    penalty <- c(penalty, rep(0, nlevel - 2), 1)
  }
  return(c(0,penalty))
}

# Main function to set up and fit the penalized GLM model - non-interaction

#' strat_poisson
#'
#' @description Fits the LASSO tree model on the given data
#' @param df data frame with variables of interest
#' @param y response variable (poisson distributed)
#' @param vars character vector of variable names to include in the model
#' @param lambda tuning parameter for the LASSO penalty
#' @param tol tolerance level for convergence
#' @param max_iter maximum number of iterations
#' @param warm.start initial values for beta
#' @return a vector of model coefficients
#' @examples
#' data("data")
#' strat_poisson(data, 'hospital_days', c('age_group', 'pulm_impair'), lambda = 0)
#' @export
strat_poisson <- function(df, y, vars, lambda=0, tol = 1e-6, max_iter = 100, warm.start) {
  # Convert specified columns to factors
  nlevels <- c()
  for(var in vars) {
    df[[var]] <- as.factor(df[[var]])
    nlevels <- c(nlevels, length(levels(df[[var]])))
  }

  #Create design matrix for all vars - no intercept
  formula_str <- paste("~", paste(vars, collapse = " + "))

  X <- model.matrix(as.formula(formula_str), data = df)

  # Total number of beta coefficients to estimate
  length_beta <- ncol(X)  # Use total columns in X for testing

  # Create the constraint matrix
  Amat <- get_A(nlevels)
  bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0

  y <- df[[y]]
  log_Y <- log(y + 0.1)

  # Initial coefficient estimates

  if (missing(warm.start)) {
    beta <- rep(1, ncol(X))
  } else {
    beta <- warm.start
  }

  #beta <- rep(1, ncol(X))
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
    solution <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec)

    beta_new_star <- solution$solution  # Extract solution

    # Check convergence
    epsilon <- sqrt(t(beta_new_star - beta) %*% (beta_new_star - beta))
    beta <- beta_new_star  # Update beta

    #iterate
    iter <- iter + 1
  }

  return(beta)
}

#Create function to find BIC for a given model

#' bic_poisson_simp
#'
#' @description Find the BIC for a given LASSO tree interaction model
#' @param df data frame with variables of interest
#' @param y response variable (poisson distributed)
#' @param vars character vector of variable names to include in the model
#' @param beta vector of model coefficients
#' @return BIC value
#' @examples
#' data("data")
#' bic_poisson_simp(data, 'hospital_days', c('age_group', 'pulm_impair'), beta = rep(1, 10))
#' @export
bic_poisson_simp <- function(df, y, vars, beta) {

  for(var in vars) {
    df[[var]] <- as.factor(df[[var]])
  }

  #Create design matrix for all vars
  formula_str <- paste("~", paste(vars, collapse = " + "))

  X <- model.matrix(as.formula(formula_str), data = df)
  y <- df[[y]]

  n <- nrow(df)
  eta <- X %*% beta
  mu <- exp(eta)
  loglik <- sum(y * eta - mu - lfactorial(y))
  bic <- -2 * loglik + log(n) * length(unique(round(beta,8)))
  return(bic)
}


#Find optimal lambda over grid using BIC

#' bic_mse_lambda_simp
#' @description Find the optimal lambda value over a grid using BIC
#' @param df data frame with variables of interest
#' @param y response variable (poisson distributed)
#' @param vars character vector of variable names to include in the model
#' @param lambda_grid grid of lambda values to test
#' @return a list containing the optimal lambda value, BIC values, and MSE values
#' @examples
#' data("data")
#' bic_mse_lambda_simp(data, 'hospital_days', vars = c('age_group', 'pulm_impair'), seq(0.1, 1, 0.1))
#' @export
bic_mse_lambda_simp <- function(df, y, vars, lambda_grid){
  bic_values <- c()
  n_unique_beta <- c()
  MSE <- c()

  for(var in vars) {
    df[[var]] <- as.factor(df[[var]])
  }

  y1 <- df[[y]]

  #Create design matrix for all vars
  formula_str <- paste("~", paste(vars, collapse = " + "))

  X <- model.matrix(as.formula(formula_str), data = df)


  for(i in 1:length(lambda_grid)){
    if(i == 1){
      beta <- strat_poisson(df, y, vars, lambda = lambda_grid[i], tol = 1e-8, max_iter = 100)
    }else{
      beta <- strat_poisson(df, y, vars, lambda = lambda_grid[i], tol = 1e-8, max_iter = 100, warm.start = beta)
    }

    bic_values <- c(bic_values, bic_poisson_simp(df, y, vars, beta))
    n_unique_beta <- c(n_unique_beta, length(unique(beta)))
    MSE <- c(MSE, mean((exp(X %*% beta) - y1)^2))
  }
  result <- list(lambda.min = lambda_grid[which.min(bic_values)], bic_values = bic_values, MSE = MSE)

  return(result)
}



## Interaction model functions

#Constraint matrix function

#' get_Ax
#'
#' @description Creates the matrix of constraints to be passed into quadratic programming problem for interaction model, code adapted from the glidars R package
#' @param p number of levels for the first categorical variable
#' @param q number of levels for the second categorical variable
#' @return An n x p matrix of constraints, where n is the number of constraints, and p is the length of beta
#' @examples
#' get_Ax(7, 4)
#' @export
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
  Ax <- rbind(c(1, rep(0, ncol(Ax)-1)), Ax)
  return(Ax)
}

#Penalty vector function

#' penalty_vec
#'
#' @description Creates the matrix of constraints to be passed into quadratic programming problem for the interaction model
#' @param p number of levels for the first categorical variable
#' @param q number of levels for the second categorical variable
#' @return A vector of length beta that when multiplied by beta gives the LASSO penalty
#' @examples
#' penalty_vec(7, 4)
#' @export
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

#' strat_poissonX
#' @description Fits the LASSO tree model on the given data assuming an interaction grid
#' @param df data frame with variables of interest
#' @param y response variable (poisson distributed)
#' @param var1 first categorical variable
#' @param var2 second categorical variable
#' @param lambda tuning parameter for the LASSO penalty
#' @param tol tolerance level for convergence
#' @param max_iter maximum number of iterations
#' @param warm.start initial values for beta
#' @return a vector containing the beta coefficients for the interaction grid
#' @examples
#' data("data")
#' strat_poissonX(data, 'hospital_days', var1 = 'age_group', var2 = 'pulm_impair', lambda = 0)
#' @export
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
    solution <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec, meq = 0)

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


#Create function to find BIC for a given model

#' bic_poisson
#'
#' @description Find the BIC for a given LASSO tree interaction model
#' @param df data frame with variables of interest
#' @param y response variable (poisson distributed)
#' @param var1 first categorical variable
#' @param var2 second categorical variable
#' @param beta vector of model coefficients
#' @return BIC value
#' @examples
#' data("data")
#' bic_poisson(data, 'hospital_days', var1 = 'age_group', var2 = 'pulm_impair', beta = rep(1, 28))
#' @export
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

#' bic_mse_lambda
#' @description Find the optimal lambda value over a grid using BIC
#' @param df data frame with variables of interest
#' @param y response variable (poisson distributed)
#' @param var1 first categorical variable
#' @param var2 second categorical variable
#' @param lambda_grid grid of lambda values to test
#' @return a list containing the optimal lambda value, BIC values, and MSE values
#' @examples
#' data("data")
#' bic_mse_lambda(data, 'hospital_days', var1 = 'age_group', var2 = 'pulm_impair', seq(0.1, 1, 0.1))
#' @export
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


