#' Simulates a matrix of genotypes for a single hypothetical variant
#'
#' @param N Number of individuals
#' @param K Number of measurement timepoints (assumes equal number per individual for now)
#' @param maf Genotype minor allele frequency
#'
#' @returns Genotype matrix of dimension N_individuals x K timepoints
#'
#' @export

simulate_genotype <- function(N, K, maf=0.25) {
  G_vec <- rbinom(N, 2, maf)
  G_mat <- matrix(G_vec, nrow=N, ncol=K, byrow=FALSE)
  G_mat
}


#' Simulates a matrix of time-varying exposures
#'
#' @param N Number of individuals
#' @param icc_E Intraclass correlation coefficient for the exposure (clustering within an individual across time)
#' @param var_e_E Error variance for the exposure
#' @param K Number of measurement timepoints (assumes equal number per individual for now)
#'
#' @returns Exposure matrix of dimension N_individuals x K timepoints
#'
#' @export

simulate_exposure <- function(N, K, icc_E, var_e_E) {
  cov_E <- var_e_E * icc_E
  sigma_E <- matrix(cov_E, nrow=K, ncol=K)
  diag(sigma_E) <- var_e_E
  E_mat <- MASS::mvrnorm(N, rep(0, K), sigma_E)
  E_mat
}


#' Simulates a longitudinal outcome variable
#'
#' @param t_mat Matrix of times for each measurement
#' @param G_mat Matrix of genotypes (N x K)
#' @param E_mat Matrix of exposures (N x K)
#' @param beta_tY Linear slope of time effect on Y
#' @param beta_GY Genotype effect size on Y
#' @param beta_EY Exposure effect size on Y
#' @param icc_Y Intraclass correlation coefficient for Y (clustering within an individual across time)
#' @param var_e_Y Error variance for Y
#'
#' @returns Outcome matrix of dimension N_individuals x K timepoints
#'
#' @export

simulate_outcome <- function(t_mat, G_mat, E_mat,
                             beta_tY, beta_GY, beta_EY,
                             icc_Y, var_e_Y) {

  # For now, assuming a compound symmetry covariance pattern

  K <- ncol(t_mat)
  stopifnot(K == ncol(E_mat))

  mu_Y <- t_mat * beta_tY + G_mat * beta_GY + E_mat * beta_EY

  cov_Y <- var_e_Y * icc_Y
  sigma_Y <- matrix(cov_Y, nrow=K, ncol=K)
  diag(sigma_Y) <- var_e_Y
  Y_mat <- apply(mu_Y, 1, function(mu_vec) MASS::mvrnorm(1, mu_vec, sigma_Y)) %>%
    t()
  Y_mat
}
