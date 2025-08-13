#' Posterior Summary for Model Fit
#'
#' Compute posterior summaries (mean and credible intervals) for estimated
#' probabilities and latent class assignments from a fitted model object.
#'
#' @param fit A fitted model object containing model_data and sample_list as
#' returned by the model fitting procedure.
#' @param cred_mass Numeric. The mass of the credible interval (default: 0.95).
#' @param keep_fit Logical. If TRUE, the original fit object is
#' included in the output (default: FALSE).
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{prob_summ_df}{A data frame with posterior mean and credible
#'     interval for probabilities, by time and class.}
#'     \item{model_data}{The model data object used for fitting.}
#'     \item{nw_sample}{Number of unique classes sampled at each iteration.}
#'     \item{posterior_mean}{List with posterior means for lambda and alpha.}
#'   }
#'   If keep_fit = TRUE, the original fit object is also included.
#'
#' @details
#' This function extracts posterior samples from the fitted model, calculates the
#' posterior mean and highest density intervals (HDI) for probabilities,
#' and summarizes the latent class assignments.
#' The summary is returned as a data frame along with class assignment
#' estimates and other relevant information. Posterior means for lambda and
#' alpha are also included if available.
#'
#' @importFrom dplyr bind_cols select filter
#' @importFrom HDInterval hdi
#' @export
#' @examples
#' \dontrun{
#' # Assume 'fit' is a fitted model object returned by your fitting function
#' summary <- posterior_summary(fit, cred_mass = 0.9)
#' head(summary$prob_summ_df)
#' }
posterior_summary = function(fit, cred_mass = 0.95, keep_fit = FALSE) {

  M = fit$model_data$M
  sample = fit$sample_list
  n_time = fit$model_data$n_time

  prob_mean = sample$prob |> comp_post_stat(mean)
  prob_lower = sample$prob |> comp_post_stat(stat_function = function(x){HDInterval::hdi(x)[1]})
  prob_upper = sample$prob |> comp_post_stat(stat_function = function(x){HDInterval::hdi(x)[2]})

  time = fit$model_data$id_time_df$time |> unique() |> rep(times = M)
  w = (1:M) |> rep(each = n_time)

  w_est = fit$w_class
  nw_sample = fit$sample_list$w |> apply(MARGIN = 1, FUN = function(x) length(unique(x)))

  prob_summ_df = dplyr::bind_cols(
    prob_mean |> trans_summ_prob(col_name = "mean", time = time, w = w),
    prob_lower |> trans_summ_prob(col_name = "li", time = time, w = w) |> dplyr::select(li),
    prob_upper |> trans_summ_prob(col_name = "ui", time = time, w = w) |> dplyr::select(ui)
  ) |>
    dplyr::filter(w %in% unique(w_est))

  if(length(dim(sample$lambda)) == 3) {
    lambda_mean = sample$lambda |> comp_post_stat(mean)
  }else{
    lambda_mean = NULL
  }

  alpha_mean = sample$alpha |> comp_post_stat(mean)

  out = list(
    prob_summ_df = prob_summ_df,
    model_data = fit$model_data,
    nw_sample = nw_sample,
    posterior_mean = list(
      lambda = lambda_mean,
      alpha = alpha_mean
    )
  )

  if(keep_fit == TRUE) {
    model_data$fit = fit
  }

  return(out)
}

