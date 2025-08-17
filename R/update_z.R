#' Update Primary Cluster Assignments (z)
#'
#' Internal: Update the primary cluster assignments (`z`) for each observation
#' using Gibbs sampling based on the current parameters of the observation
#' model (mu and sigma) and optionally the sequence model (w and alpha).
#'
#' @param mu Numeric matrix of cluster means (G x number of variables).
#' @param sigma Numeric vector of cluster standard deviations (length G).
#' @param w Optional integer vector of secondary cluster assignments (default NULL).
#' @param alpha Optional numeric matrix of coefficients for the sequence model (default NULL).
#' @param model_data List containing model data including `x`, `G`, `n_vars`, `id_time_df`.
#'
#' @return A list with:
#'   \describe{
#'     \item{z}{Integer vector of sampled primary cluster assignments.}
#'     \item{z_post_prob}{Matrix of posterior probabilities for each observation belonging to each primary cluster.}
#'   }
#' @importFrom stats dnorm
#' @importFrom extraDistr rcatlp
#' @importFrom mclust logsumexp
#' @keywords internal
update_z = function(mu, sigma, w = NULL, alpha = NULL, model_data) {

  x = model_data$x
  G = model_data$G
  n = nrow(x)
  n_vars = model_data$n_vars
  id = model_data$id_time_df$id
  time_seq = model_data$id_time_df$time_seq

  log_like = matrix(NA, nrow = n, ncol = G)

  for (g in 1:G) {
    log_like[, g] = rowSums(
      dnorm(x,
            mean = matrix(mu[g, ], nrow = n, ncol = n_vars, byrow = T),
            sd = sigma[g], log = TRUE)
    )

    if(!is.null(w)) {

      pzw = predict_cat(
        z = rep(g, n),
        alpha = alpha,
        w_vec = w[id],
        time_seq = time_seq,
        model_data = model_data,
        type = "category probability",
        intercept = FALSE
      )

      log_like[, g] = log_like[, g] + log(pzw)

    }
  }

  # Normalize log-likelihoods to get log-probabilities
  log_like = log_like - matrix(mclust::logsumexp(log_like), nrow = n, ncol = G, byrow = F)

  # Sample new z from the categorical
  z = extraDistr::rcatlp(n = n, log_prob = log_like) + 1

  out = list(z = z, z_post_prob = exp(log_like))
  return(out)

}




