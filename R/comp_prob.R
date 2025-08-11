#' Compute Cluster Probabilities from Linear Predictors
#'
#' Internal: Compute cluster membership probabilities using the softmax transformation of the linear predictors for each observation.
#'
#' @param alpha Numeric matrix of coefficients (basis × groups).
#' @param model_data List. Model data structure containing X_unique and n_id.
#'
#' @return Numeric matrix of probabilities for each observation and group.
#'
#' @details
#' This function computes the linear predictors for each observation using the design matrix X_unique and the coefficient matrix alpha, then applies the softmax transformation (from mclust) to obtain cluster membership probabilities for each observation.
#'
#' @importFrom mclust softmax
#' @keywords internal
comp_prob = function(alpha, model_data) {
  X_unique = model_data$X_unique
  n_id = model_data$n_id

  prob = mclust::softmax(
    X_unique %*% alpha
  )

  return(prob)
}
