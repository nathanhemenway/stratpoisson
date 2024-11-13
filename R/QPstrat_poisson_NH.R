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
    constraint <- rep(0, length_beta)
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

  # Combine all constraints
  Amat <- do.call(rbind, c(non_neg_constraints, monotonicity_constraints))

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
  beta <- c(0, rep(1, ncol(X)-1))
  epsilon <- 99
  iter <- 0

  #IRWLS
  while (epsilon > tol & iter <= max_iter) {
    eta <- X %*% beta
    mu <- exp(eta)
    nu <- exp(eta)

    A <- .sparseDiagonal(as.vector(nu), n = nrow(nu))
    z <- eta + solve(A) %*% (y - mu)

    # Compute D and d for the quadratic program
    D <- 2 * t(X) %*% A %*% X  # Quadratic term
    d <- 2 * t(X) %*% A %*% z - lambda * (beta[tail(beta_indices[[1]], 1)] + beta[tail(beta_indices[[2]], 1)])  # Linear term with penalty

    # Solve the quadratic program with constraints
    solution <- solve.QP(Dmat = D, dvec = d, Amat = t(Amat), bvec = bvec, meq = 1)

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
  Ax <- rbind(c(1, rep(0, p*q-1)), Ax)
  return(Ax)
}


strat_poissonX <- function(df, y, var1, var2, lambda=0, tol = 1e-6, max_iter = 100) {
  # Convert specified columns to factors
  df[[var1]] <- as.factor(df[[var1]])
  df[[var2]] <- as.factor(df[[var2]])

  # Find the number of levels in var1 and var2
  nlevels1 <- length(levels(df[[var1]]))
  nlevels2 <- length(levels(df[[var2]]))

  print(paste0('p = ', nlevels1, ', q = ', nlevels2))

  #Create interaction terms
  df$interaction <- interaction(df[[var1]], df[[var2]])

  # Initialize the design matrix
  X <- model.matrix(~interaction, data = df)

  # Create the constraint matrix
  Amat <- get_Ax(p = nlevels1, q = nlevels2)
  bvec <- rep(0, nrow(Amat))  # All constraints are of the form >= 0

  y <- df[[y]]

  # Initial coefficient estimates
  #beta <- solve(t(X) %*% X) %*% t(X) %*% log_Y
  beta <- c(0, rep(1, ncol(X)-1))
  epsilon <- 99
  iter <- 0

  #IRWLS
  while (epsilon > tol & iter <= max_iter) {
    eta <- X %*% beta
    mu <- exp(eta)
    nu <- exp(eta)

    A <- .sparseDiagonal(as.vector(nu), n = nrow(nu))
    z <- eta + solve(A) %*% (y - mu)

    # Compute D and d for the quadratic program
    D <- 2 * t(X) %*% A %*% X  # Quadratic term
    d <- 2 * t(X) %*% A %*% z - lambda * (beta[tail(beta_indices[[1]], 1)] + beta[tail(beta_indices[[2]], 1)])  # Linear term with penalty

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


##Test out the functions

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
  select(hospital_days, age_group, pulm_impair, age_pulm) %>%
  drop_na() -> data

resX <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = 0.1, tol = 1e-8, max_iter = 100)

matrix(resX, nrow = 7, ncol=4)

library(pheatmap)
pheatmap(matrix(resX, nrow = 7, ncol=4), cluster_rows = F, cluster_cols = F)
