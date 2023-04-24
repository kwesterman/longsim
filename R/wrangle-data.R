#' Converts longitudinal data from "wide" to "long" format
#' (assumes that the incoming dataset is stored as a matrix object)
#'
#' @param wide_mat Matrix of dimension N_individuals x K_timepoints.
#' @param var_name Variable name for the data values
#'
#' @returns Data frame containing longitudinal data in long format
#'
#' @export

convert_wide_mat_to_long_df <- function(wide_mat, var_name) {
  colnames(wide_mat) <- paste0("t", seq(1, ncol(wide_mat)))
  id_vec <- if (!is.null(rownames(wide_mat))) rownames(wide_mat) else seq(1, nrow(wide_mat))
  tibble::as_tibble(wide_mat, .name_repair="unique") %>%
    mutate(id = id_vec)  %>%
    tidyr::pivot_longer(-id, names_to="timept", values_to=var_name)
}
