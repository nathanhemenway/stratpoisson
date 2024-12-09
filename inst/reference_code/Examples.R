#Examples
#Test out the functions
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
    grid[age_group_levels == age, pulm_impair_levels == pulm] <- predict(fit, newdata = data.frame('Age.Group' = as.ordered(age),
                                                                                                   'Pulmonary.Impairment' = as.ordered(pulm)))
  }
}

#Create heatmap
pheatmap(grid, cluster_rows = F, cluster_cols = F, display_numbers = T,
         show_rownames = T, show_colnames = T, fontsize = 8,
         labels_row = age_group_levels,
         labels_col = pulm_impair_levels,
         angle_col=45,
         main = 'Estimated Average Hospitalization Days by Age Group and Pulmonary Impairment Status')


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


