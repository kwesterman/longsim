% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate-data.R
\name{simulate_exposure}
\alias{simulate_exposure}
\title{Simulates a matrix of time-varying exposures}
\usage{
simulate_exposure(
  N,
  K,
  icc_E,
  var_e_E,
  C_mat = NULL,
  beta_CE = NULL,
  frac_missing = 0
)
}
\arguments{
\item{N}{Number of individuals}

\item{K}{Number of measurement timepoints (assumes equal number per individual for now)}

\item{icc_E}{Intraclass correlation coefficient for the exposure (clustering within an individual across time)}

\item{var_e_E}{Error variance for the exposure}

\item{C_mat}{Matrix of confounders (N x K)}

\item{beta_CE}{Confounder effect size on E}

\item{frac_missing}{Fraction of all values set to missing (MCAR)}
}
\value{
Exposure matrix of dimension N_individuals x K timepoints
}
\description{
Simulates a matrix of time-varying exposures
}
