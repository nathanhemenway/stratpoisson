#test code here

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
