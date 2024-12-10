## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(stratpoisson)

## -----------------------------------------------------------------------------
res1 <- strat_poisson(data, 'hospital_days', c('age_group', 'pulm_impair', 'num_comorb'), lambda = 0, tol = 1e-8, max_iter = 100)
res1

## -----------------------------------------------------------------------------
lam_result1 <- bic_mse_lambda_simp(df=data, y = 'hospital_days', vars = c('age_group', 'pulm_impair', 'num_comorb'), lambda_grid=seq(from=0, to=5, by=0.1))
plot(x = seq(from=0, to=5, by=0.1), y = lam_result1$bic_values, type='l')
lambda1 <- lam_result1$lambda.min
lambda1 #For this data the best lambda = 0

## -----------------------------------------------------------------------------
resX <- strat_poissonX(data, 'hospital_days', 'age_group', 'pulm_impair', lambda = 0, tol = 1e-8, max_iter = 100)

## -----------------------------------------------------------------------------
matrix(exp(resX), nrow = length(levels(data$age_group)))

## -----------------------------------------------------------------------------
library(pheatmap)
pheatmap(matrix(exp(resX), nrow = length(levels(data$age_group))), cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = c('<= 29', '30-39', '40-49', '50-59', '60-69', '70-79', '80+'),
         labels_col = c('None', 'Mild', 'Moderate', 'Severe'),
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status')

## -----------------------------------------------------------------------------
lam_result <- bic_mse_lambda(df=data, y = 'hospital_days', var1 = 'age_group', var2 = 'pulm_impair', lambda_grid=seq(from=0, to=5, by=0.1))
plot(x = seq(from=0, to=5, by=0.1), y = lam_result$bic_values, type='l')
lambda1 <- lam_result$lambda.min
lambda1 #For this data the best lambda = 0

