#' Simulates a matrix of genotypes for a single hypothetical variant
#'
#' @param N Number of individuals
#' @param K Number of measurement timepoints (assumes equal number per individual for now)
#' @param maf Genotype minor allele frequency
#' @param frac_missing Fraction of all values set to missing (MCAR)
#'
#' @returns Genotype matrix of dimension N_individuals x K timepoints
#'
#' @export

simulate_genotype <- function(N, K, maf=0.25, frac_missing=0) {
  G_vec <- rbinom(N, 2, maf)
  G_mat <- matrix(G_vec, nrow=N, ncol=K, byrow=FALSE)

  missing_idx <- sample(c(TRUE, FALSE), prob=c(frac_missing, 1 - frac_missing),
                        size=length(G_mat), replace=TRUE)
  G_mat[missing_idx] <- NA

  G_mat
}


#' Simulates a matrix of time-varying covariates
#'
#' @param N Number of individuals
#' @param K Number of measurement timepoints (assumes equal number per individual for now)
#' @param icc_C Intraclass correlation coefficient for the covariate (clustering within an individual across time)
#' @param var_e_C Error variance for the covariate
#' @param frac_missing Fraction of all values set to missing (MCAR)
#'
#' @returns Covariate matrix of dimension N_individuals x K timepoints
#'
#' @export

simulate_covariate <- function(N, K, icc_C, var_e_C, frac_missing=0) {
  cov_C <- var_e_C * icc_C
  sigma_C <- matrix(cov_C, nrow=K, ncol=K)
  diag(sigma_C) <- var_e_C
  C_mat <- MASS::mvrnorm(N, rep(0, K), sigma_C)

  missing_idx <- sample(c(TRUE, FALSE), prob=c(frac_missing, 1 - frac_missing),
                        size=length(C_mat), replace=TRUE)
  C_mat[missing_idx] <- NA

  C_mat
}


#' Simulates a matrix of time-varying exposures
#'
#' @param N Number of individuals
#' @param K Number of measurement timepoints (assumes equal number per individual for now)
#' @param icc_E Intraclass correlation coefficient for the exposure (clustering within an individual across time)
#' @param var_e_E Error variance for the exposure
#' @param C_mat Matrix of confounders (N x K)
#' @param beta_CE Confounder effect size on E
#' @param frac_missing Fraction of all values set to missing (MCAR)
#'
#' @returns Exposure matrix of dimension N_individuals x K timepoints
#'
#' @export

simulate_exposure <- function(N, K, icc_E, var_e_E, C_mat=NULL, beta_CE=NULL, frac_missing=0) {
  cov_E <- var_e_E * icc_E
  sigma_E <- matrix(cov_E, nrow=K, ncol=K)
  diag(sigma_E) <- var_e_E

  if (is.null(C_mat)) {
    E_mat <- MASS::mvrnorm(N, rep(0, K), sigma_E)
  } else {
    mu_E <- C_mat * beta_CE
    E_mat <- apply(mu_E, 1, function(mu_vec) MASS::mvrnorm(1, mu_vec, sigma_E)) %>%
      t()
  }

  missing_idx <- sample(c(TRUE, FALSE), prob=c(frac_missing, 1 - frac_missing),
                        size=length(E_mat), replace=TRUE)
  E_mat[missing_idx] <- NA

  E_mat
}


#' Simulates a longitudinal outcome variable
#'
#' @param matrix_list List of matrices affecting the mean of Y
#' @param beta_list List of effect sizes associated with the matrices in matrix_list
#' @param icc_Y Intraclass correlation coefficient for Y (clustering within an individual across time)
#' @param var_e_Y Error variance for Y
#' @param frac_missing Fraction of all values set to missing (MCAR)
#'
#' @returns Outcome matrix of dimension N_individuals x K timepoints
#'
#' @export

simulate_outcome <- function(matrix_list, beta_list, icc_Y, var_e_Y, frac_missing=0) {

  # For now, assuming a compound symmetry covariance pattern

  N <- nrow(matrix_list[[1]])
  K <- ncol(matrix_list[[1]])

  mu_Y <- matrix(0, nrow=N, ncol=K)
  stopifnot(length(matrix_list) == length(beta_list))
  for (i in seq(1, length(matrix_list))) {
    mu_Y <- mu_Y + matrix_list[[i]] * beta_list[[i]]
  }

  cov_Y <- var_e_Y * icc_Y
  sigma_Y <- matrix(cov_Y, nrow=K, ncol=K)
  diag(sigma_Y) <- var_e_Y
  Y_mat <- apply(mu_Y, 1, function(mu_vec) MASS::mvrnorm(1, mu_vec, sigma_Y)) %>%
    t()

  missing_idx <- sample(c(TRUE, FALSE), prob=c(frac_missing, 1 - frac_missing),
                        size=length(Y_mat), replace=TRUE)
  Y_mat[missing_idx] <- NA

  Y_mat
}


#' Simulate a longitudinal dataset
#'
#' @param N Number of individuals
#' @param K Number of measurement timepoints (assumes equal number per individual for now)
#'
#' @param maf Genotype minor allele frequency
#'
#' @param icc_C Intraclass correlation coefficient for the covariate (clustering within an individual across time)
#' @param var_e_C Error variance for the covariate
#'
#' @param beta_CE Confounder effect size on E
#' @param icc_E Intraclass correlation coefficient for the exposure (clustering within an individual across time)
#' @param var_e_E Error variance for the exposure
#'
#' @param beta_tY Time effect on Y
#' @param beta_GY Genotype effect on Y
#' @param beta_CY Covariate effect on Y
#' @param beta_EY Exposure effect on Y
#' @param icc_Y Intraclass correlation coefficient for Y (clustering within an individual across time)
#' @param var_e_Y Error variance for Y
#'
#' @param frac_missing Fraction of all values set to missing (MCAR)
#'
#' @returns Data frame in long format containing fields for ID (id), timepoint label (timept),
#' outcome (Y), measurement time (t), genotype (G),exposure (E), and covariate (C)
#'
#' @export

simulate_long_data <- function(
    N = 100,
    K = 3,
    maf = 0.25,
    icc_C = 0.5,
    var_e_C = 1,
    beta_CE = 0,
    icc_E = 0.5,
    var_e_E = 1,
    beta_tY = -1,
    beta_GY = 0.5,
    beta_CY = 0,
    beta_EY = 1,
    icc_Y = 0.8,
    var_e_Y = 1,
    frac_missing = 0
) {
  tmat <- matrix(seq(1, K), nrow=N, ncol=K, byrow=TRUE)
  gmat <- simulate_genotype(N, K, maf=maf, frac_missing=frac_missing)
  cmat <- simulate_covariate(N, K, icc_C, var_e_C, frac_missing=frac_missing)
  emat <- simulate_exposure(N, K,
                            C_mat=cmat, beta_CE=beta_CE,
                            icc_E, var_e_E,
                            frac_missing=frac_missing)
  ymat <- simulate_outcome(list(tmat, gmat, cmat, emat),
                           list(beta_tY, beta_GY, beta_CY, beta_EY),
                           icc_Y, var_e_Y, frac_missing=frac_missing)

  merge_variables(list(Y = ymat, t = tmat, G = gmat, E = emat, C = cmat)) %>%
    mutate(id = factor(id),
           timept = factor(timept, levels=paste0("t", seq(1, K))))
}
