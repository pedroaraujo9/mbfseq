#' Update Group-Specific Coefficients via Polya-Gamma Augmentation
#'
#' Internal: Compute the conditional posterior mean and draw a sample of the group-specific coefficient vector alpha for the g-th group in a multinomial logistic spline model, using Polya-Gamma data augmentation.
#'
#' @param g Integer. Index of the group currently being updated.
#' @param Z Matrix. Binary indicator matrix for cluster membership (rows = observations, columns = clusters).
#' @param X Matrix. Design matrix for fixed and random effects.
#' @param linear_pred Matrix. Current linear predictors for all groups.
#' @param inv_cov_list List. List of inverse prior covariance matrices for each group.
#'
#' @return Numeric vector. A posterior draw of the coefficient vector for group g.
#'
#' @details
#' For the specified group g, this function computes the Polya-Gamma augmented weights, constructs the conditional posterior for the group-specific coefficients, and draws a sample using the sample_alpha helper. Used internally in the Gibbs sampler for multinomial logistic models with basis expansions.
#'
#' @importFrom BayesLogit rpg
#' @keywords internal
augment_alpha = function(g,
                         Z,
                         X,
                         linear_pred,
                         inv_cov_list) {

  max_other = apply(linear_pred[, -g, drop = FALSE], 1, max)
  C = max_other + log(rowSums(exp(linear_pred[, -g, drop = FALSE] - max_other)))
  n = nrow(Z)

  omega = BayesLogit::rpg(n, h = 1, z = linear_pred[, g, drop = FALSE] - C)
  center = matrix(0, nrow = ncol(X))

  alpha = sample_alpha(
    X = X,
    omega = omega,
    inv_cov = inv_cov_list[[g]],
    z = Z[, g, drop = FALSE],
    C = C,
    center = center
  ) |> as.numeric()

  return(alpha)
}




