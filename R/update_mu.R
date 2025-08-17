#' Update Cluster Means (mu)
#'
#' Internal: Update the cluster means (`mu`) using Gibbs sampling for the observation model.
#'
#' @param z Integer vector. The current cluster assignments for each observation.
#' @param sigma Numeric vector. The current cluster standard deviations.
#' @param model_data List. Model data structure containing `G`, `x`, and `n_vars`.
#' @param mu_fixed_sd Numeric. Standard deviation for the prior on the cluster means (default: 1).
#'
#' @return Numeric matrix. A matrix where each row represents the updated mean for a cluster and columns represent variables.
#'
#' @details
#' This function updates the mean for each cluster based on the observations assigned to that cluster (`x[z == g, ]`), the current cluster standard deviation, and a fixed prior standard deviation (`mu_fixed_sd`). The update is performed independently for each cluster and each variable using a normal distribution derived from the conditional posterior.
#'
#' @importFrom stats rnorm
#' @keywords internal
update_mu = function(z, sigma, model_data, mu_fixed_sd = 1) {

  G = model_data$G
  x = model_data$x
  n_vars = model_data$n_vars

  mu = matrix(NA, nrow = G, ncol = n_vars)
  ng = z %>% factor(levels = 1:G) %>% table()

  for (g in 1:G) {
    var_g = 1/((1/(mu_fixed_sd^2)) + ng[g]/(sigma[g]^2))
    mu_g = var_g * (ng[g] * colMeans(x[z == g, ])/(sigma[g]^2))
    mu[g, ] = rnorm(n_vars, mean = mu_g, sd = sqrt(var_g))
  }

  return(mu)

}
