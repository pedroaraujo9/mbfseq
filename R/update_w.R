#' Update Latent Cluster Memberships
#'
#' Internal: Sample updated latent cluster assignments for each individual based on posterior probabilities, given current model parameters and observed data.
#'
#' @param alpha Matrix. Current group-specific coefficient matrix (columns = groups).
#' @param pw Numeric vector. Prior cluster membership probabilities for each cluster.
#' @param z Integer vector. Observed cluster labels for each observation.
#' @param model_data List. Model data structure (from create_model_data), must include M, n_id, id_unique, w_logp_base.
#' @param intercept Logical. Include intercept in model? Default is FALSE.
#' @param mean_effect Logical. Use mean effect? Default is FALSE.
#'
#' @return A list with:
#'   \describe{
#'     \item{w}{Integer vector of sampled cluster assignments.}
#'     \item{w_post_prob}{Matrix of posterior cluster probabilities.}
#'     \item{marg_log_like}{Vector of marginal log-likelihood values for each individual.}
#'   }
#'
#' @details
#' For each individual, this function computes the posterior distribution over latent cluster indicators using the current model parameters. It:
#' \itemize{
#'   \item Computes the log-probability of the observed data for each possible cluster assignment.
#'   \item Adds the log prior for each cluster.
#'   \item Normalizes to obtain posterior probabilities (in log space for stability).
#'   \item Samples new cluster assignments from the posterior using categorical sampling.
#'   \item Returns the sampled assignments, the posterior probability matrix, and the marginal log-likelihood for each individual.
#' }
#'
#' @importFrom dplyr group_by summarise mutate ungroup select left_join
#' @importFrom tidyr spread
#' @importFrom extraDistr rcatlp
#' @keywords internal
update_w = function(alpha,
                    pw,
                    z,
                    model_data,
                    intercept = FALSE,
                    mean_effect = FALSE) {

  M = model_data$M
  id_unique = model_data$id_unique
  n_id = model_data$n_id
  w_logp_base = model_data$w_logp_base

  w_logp_base$z = rep(z, times = M)

  w_logp_base$logp_z = log(predict_cat(
    alpha = alpha,
    z = w_logp_base$z,
    w_vec = w_logp_base$w,
    time_seq = w_logp_base$time_seq,
    model_data = model_data,
    type = "category probability",
    intercept = intercept
  ))

  logpw = log(rep(pw, times = n_id))

  sum_logp_cond_z = w_logp_base |>
    dplyr::group_by(id, w) |>
    dplyr::summarise(sum_logp_z = sum(logp_z), .groups = "drop") |>
    dplyr::ungroup()

  sum_logp_z = sum_logp_cond_z |>
    dplyr::mutate(sum_logp_z = sum_logp_z + logpw)

  sum_logp_z$sum_logp_z[sum_logp_z$sum_logp_z < -1e300] = -1e300

  marg_log_like = sum_logp_z |>
    dplyr::group_by(id) |>
    dplyr::summarise(logp_z = logsumexp_cpp(sum_logp_z), .groups = "drop") |>
    ungroup()

  logp_w = sum_logp_z |>
    dplyr::left_join(marg_log_like, by = "id") |>
    dplyr::mutate(logp_w = sum_logp_z - logp_z) |>
    dplyr::select(id, w, logp_w) |>
    tidyr::spread(w, logp_w) |>
    dplyr::select(-id) |>
    as.matrix()

  w = extraDistr::rcatlp(n = 1:n_id, log_prob = logp_w) + 1
  names(w) = id_unique

  out = list(
    w = w,
    w_post_prob = exp(logp_w),
    marg_log_like = marg_log_like$logp_z
  )

  return(out)

}



