#' @title Update Lambda Parameters
#' @description
#' Internal function to update the lambda (penalty) parameters for each group and class using a Gamma prior, based on the current values of the alpha coefficients.
#'
#' @param alpha Numeric matrix of current alpha coefficients (basis coefficients), with rows corresponding to basis functions and columns to groups/classes.
#' @param a Numeric. Shape parameter for the Gamma prior (default: 1).
#' @param b Numeric. Rate parameter for the Gamma prior (default: 1).
#' @param model_data A \code{model_data} object as returned by \code{\link{create_model_data}}, containing \code{G}, \code{M}, \code{n_basis}, and \code{notpen_index}.
#'
#' @details
#' For each class (except the last), this function aggregates the squared alpha coefficients (excluding non-penalized indices) by group, and samples new lambda values from the inverse Gamma posterior. The update is performed for each group and class.
#'
#' @return
#' A matrix of updated lambda values with dimensions \code{M} (groups) by \code{G-1} (classes).
#'
#' @importFrom magrittr %>%
#' @importFrom stats aggregate rgamma
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Example with dummy data
#' model_data <- list(G = 3, M = 2, n_basis = 5, notpen_index = c(1,2,6,7))
#' alpha <- matrix(rnorm(5*2*3), nrow = 10, ncol = 3)
#' update_lambda(alpha, a = 1, b = 1, model_data = model_data)
#' }
update_lambda_id = function(alpha, a = 1, b = 1, model_data) {

  G = model_data$G
  M = model_data$M
  n_id = model_data$n_id

  sq_sum = colSums(alpha[1:n_id, -G]^2)
  n_params = n_id

  new_lambda = 1/stats::rgamma(n = G-1, a + n_params/2, b + sq_sum/2)

  return(new_lambda)

}
