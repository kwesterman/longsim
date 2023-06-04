#' Converts longitudinal data from "wide" to "long" format
#' (assumes that the incoming dataset is stored as a matrix object)
#'
#' @param wide_mat Matrix of dimension N_individuals x K_timepoints
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


#' Converts a list of matrices into a single dataset in long format
#'
#' @param mat_list List of matrices, all of dimension N_individuals x K_timepoints
#'
#' @returns Data frame in long format containing a column for each variable in
#' the input matrix list
#'
#' @export

merge_variables <- function(mat_list) {
  if (is.null(names(mat_list))) names(mat_list) <- paste0("V", seq(1, length(mat_list)))
  df_list <- purrr::imap(mat_list, ~ convert_wide_mat_to_long_df(.x, .y))
  purrr::reduce(df_list, function(x, y) dplyr::inner_join(x, y, by=c("id", "timept")))
}
