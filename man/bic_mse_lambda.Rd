% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QPstrat_poisson_NH.R
\name{bic_mse_lambda}
\alias{bic_mse_lambda}
\title{bic_mse_lambda}
\usage{
bic_mse_lambda(df, y, var1, var2, lambda_grid)
}
\arguments{
\item{df}{data frame with variables of interest}

\item{y}{response variable (poisson distributed)}

\item{var1}{first categorical variable}

\item{var2}{second categorical variable}

\item{lambda_grid}{grid of lambda values to test}
}
\value{
a list containing the optimal lambda value, BIC values, and MSE values
}
\description{
Find the optimal lambda value over a grid using BIC
}
\examples{
data("data")
bic_mse_lambda(data, 'hospital_days', var1 = 'age_group', var2 = 'pulm_impair', seq(0.1, 1, 0.1))
}
