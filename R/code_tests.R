#test code here

#Create function to return the A matrix:

get_A <- function(beta_indices){
  #beta_indices: list of indices for each group of coefficients

  # Initialize constraint matrix components
  non_neg_constraints <- list()
  monotonicity_constraints <- list()

  # Example setup for grouped monotonicity constraints
  for (group in beta_indices) {
    n_levels <- length(group)

    # 1. Non-negativity constraint for the first level of this group
    constraint <- rep(0, length(beta))
    constraint[group[1]] <- 1
    non_neg_constraints <- append(non_neg_constraints, list(constraint))

    # 2. Monotonicity constraints within the group
    for (k in 2:n_levels) {
      constraint <- rep(0, length(beta))
      constraint[group[k - 1]] <- -1
      constraint[group[k]] <- 1
      monotonicity_constraints <- append(monotonicity_constraints, list(constraint))
    }
  }

  # Combine non-negativity and monotonicity constraints
  Amat <- do.call(rbind, c(non_neg_constraints, monotonicity_constraints))
  #bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0

  return(Amat)
}

#create list of indices for each group of coefficients as an example:
beta_indices <- list(c(1, 2, 3), c(4, 5, 6, 7), c(8, 9))
beta <- unlist(beta_indices)

# Initialize constraint matrix components
non_neg_constraints <- list()
monotonicity_constraints <- list()

# Example setup for grouped monotonicity constraints
for (group in beta_indices) {
  n_levels <- length(group)

  # 1. Non-negativity constraint for the first level of this group
  constraint <- rep(0, length(beta))
  constraint[group[1]] <- 1
  non_neg_constraints <- append(non_neg_constraints, list(constraint))

  # 2. Monotonicity constraints within the group
  for (k in 2:n_levels) {
    constraint <- rep(0, length(beta))
    constraint[group[k - 1]] <- -1
    constraint[group[k]] <- 1
    monotonicity_constraints <- append(monotonicity_constraints, list(constraint))
  }
}

# Combine non-negativity and monotonicity constraints
Amat <- do.call(rbind, c(non_neg_constraints, monotonicity_constraints))
bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0


strat_poisson <- function(df, y, var1, var2) {
  # Convert specified columns to factors
  df[[var1]] <- as.factor(df[[var1]])
  df[[var2]] <- as.factor(df[[var2]])

  # Find the number of levels in var1 and var2
  nlevels1 <- length(levels(df[[var1]]))
  nlevels2 <- length(levels(df[[var2]]))


  # Initialize the design matrix
  formula_str <- paste("~", var1, "+", var2)
  X <- model.matrix(as.formula(formula_str), data = df)

  # Create beta indices vector
  beta_indices <- list(c(1), rep(1, nlevels1 - 1), rep(2, nlevels2 - 1))

  # Create the constraint matrix
  Amat <- get_A(beta_indices)

  return(Amat)
}




library(NHANES)
data("NHANES")
is.ordered(NHANES$BMI_WHO)
NHANES$BMI_WHO <- factor(NHANES$BMI_WHO, ordered = TRUE)
NHANES$BMI_WHO



#Chat GPT update
# Function to build the constraint matrix with diagnostic print statements
get_A <- function(beta_indices, length_beta) {
  # Initialize constraint matrix components
  non_neg_constraints <- list()
  monotonicity_constraints <- list()

  # Set up grouped monotonicity constraints
  for (group in beta_indices) {
    n_levels <- length(group)

    # 1. Non-negativity constraint for the first level of this group
    constraint <- rep(0, length_beta)  # Length matches total number of betas
    constraint[group[1]] <- 1
    non_neg_constraints <- append(non_neg_constraints, list(constraint))

    # 2. Monotonicity constraints within the group
    for (k in 2:n_levels) {
      constraint <- rep(0, length_beta)
      constraint[group[k - 1]] <- -1
      constraint[group[k]] <- 1
      monotonicity_constraints <- append(monotonicity_constraints, list(constraint))
    }
  }

  # Combine non-negativity and monotonicity constraints
  Amat <- do.call(rbind, c(non_neg_constraints, monotonicity_constraints))

  # Print diagnostics
  # print("beta_indices:")
  # print(beta_indices)
  # print("length_beta:")
  # print(length_beta)
  # print("Amat dimensions:")
  # print(dim(Amat))

  return(Amat)
}

# Main function to set up and fit the penalized GLM model
strat_poisson <- function(df, y, var1, var2, lambda=0, tol = 1e-6, max_iter = 100) {
  # Convert specified columns to factors
  df[[var1]] <- as.factor(df[[var1]])
  df[[var2]] <- as.factor(df[[var2]])

  # Find the number of levels in var1 and var2
  nlevels1 <- length(levels(df[[var1]]))
  nlevels2 <- length(levels(df[[var2]]))

  # Initialize the design matrix, excluding the intercept
  formula_str <- paste("~", var1, "+", var2)
  X <- model.matrix(as.formula(formula_str), data = df)

  # Create beta_indices vector: identifies which columns in X correspond to each non-reference level
  var1_levels <- seq(2, nlevels1)  # Exclude reference category
  var2_levels <- seq(nlevels1 + 1, nlevels1 + nlevels2 - 1)  # Adjusted for second factor

  beta_indices <- list(var1_levels, var2_levels)

  # Total number of beta coefficients to estimate
  length_beta <- ncol(X)  # Use total columns in X for testing

  # Create the constraint matrix
  Amat <- get_A(beta_indices, length_beta)
  bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0

  y <- df[[y]]
  log_Y <- log(y + 0.1)

  # Initial coefficient estimates
  #beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y
  beta <- rep(0, ncol(X))
  epsilon <- 99
  iter <- 0

    #IRWLS
    while (epsilon > tol & iter <= max_iter) {
      eta <- X %*% beta
      mu <- exp(eta)
      nu <- exp(eta)

      V <- diag(as.vector(nu))
      Z <- eta + solve(V) %*% (y - mu)

      # Compute D and d for the quadratic program
      D <- t(X) %*% V %*% X  # Quadratic term
      d <- t(X) %*% V %*% Z - lambda * rep(1, q)  # Linear term with penalty

      # Solve the quadratic program with constraints
      solution <- solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec, meq = 0)

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

strat_poisson(df, 'DaysMentHlthBad', 'BMI_WHO', 'AgeDecade', lambda = 0.1, tol = 1e-3, max_iter = 10)

library(tidyverse)
NHANES %>%
  select(SleepHrsNight, BMI_WHO, AgeDecade) %>%
  drop_na() -> df

var1 <- 'BMI_WHO'
var2 <- 'AgeDecade'
df[[var1]] <- as.factor(df[[var1]])
df[[var2]] <- as.factor(df[[var2]])

# Find the number of levels in var1 and var2
nlevels1 <- length(levels(df[[var1]]))
nlevels2 <- length(levels(df[[var2]]))

# Initialize the design matrix, excluding the intercept
formula_str <- paste("~", var1, "+", var2)
X <- model.matrix(as.formula(formula_str), data = df)

# Create beta_indices vector: identifies which columns in X correspond to each non-reference level
var1_levels <- seq(2, nlevels1)  # Exclude reference category
var2_levels <- seq(nlevels1 + 1, nlevels1 + nlevels2 - 1)  # Adjusted for second factor

beta_indices <- list(var1_levels, var2_levels)

# Total number of beta coefficients to estimate
length_beta <- ncol(X)  # Use total columns in X for testing

# Create the constraint matrix
Amat <- get_A(beta_indices, length_beta)
bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0

y <- df[['SleepHrsNight']]
log_Y <- log(y + 1)

# Initial coefficient estimates
beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y
dim(X)

print(c(dim(X), dim(log_Y)))
