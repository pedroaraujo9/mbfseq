#' Update Cluster Standard Deviations
#'
#' Internal: Update the cluster standard deviations (`sigma`) using Gibbs sampling.
#'
#' @param mu Matrix of cluster means (G x n_vars).
#' @param z Vector of cluster assignments for each observation (length n).
#' @param sigma_a Numeric. Shape parameter for the inverse gamma prior on sigma^2 (default 1).
#' @param sigma_b Numeric. Rate parameter for the inverse gamma prior on sigma^2 (default 1).
#' @param model_data List. Model data structure including `x` (observation matrix), `G` (number of clusters), and `n_vars` (number of variables).
#'
#' @details
#' For each cluster g, this function computes the sum of squared errors (SSE)
#' between the observations assigned to that cluster and the cluster's current
#' mean. This SSE is used along with the prior parameters (`sigma_a` and `sigma_b`)
#' to define the parameters of the inverse-gamma posterior distribution for the
#' variance (sigma^2) of that cluster. A new variance is sampled from this posterior,
#' and the standard deviation (`sigma`) is obtained by taking the square root of
#' the sampled variance.
#'
#' @importFrom stats rgamma
#'
#' @return Numeric vector of updated cluster standard deviations (length G).
#'
#' @keywords internal
update_sigma = function(mu, z, sigma_a = 1, sigma_b = 1, model_data) {
  x = model_data$x
  G = model_data$G
  n_vars = model_data$n_vars

  sigma2 = numeric(G)

  for (g in 1:G) {
    idx = which(z == g)
    n_k = length(idx)

    x_k = x[idx, , drop = FALSE]
    mu_k = matrix(mu[g, ], nrow = n_k, ncol = n_vars, byrow = TRUE)

    sse = sum((x_k - mu_k)^2)

    shape = sigma_a + (n_k * n_vars) / 2
    rate = sigma_b + sse / 2

    sigma2[g] = 1 / rgamma(1, shape = shape, rate = rate)
  }

  sigma = sqrt(sigma2)
  return(sigma)
}
