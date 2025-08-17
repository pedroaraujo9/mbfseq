 #' Posterior Summary for Cluster Sequence Model
 #'
 #' Summarizes posterior samples for cluster sequence models, including means and credible intervals for stage probabilities.
 #'
 #' @param fit Output from fit_mbfseq containing posterior samples and model metadata.
 #' @param M Optional. Number of mixture components. If NULL, extracted from fit.
 #' @param G Optional. Number of clusters. If NULL, extracted from fit.
 #' @param cred_mass Credible mass for intervals (default: 0.95).
 #' @param keep_fit Logical. If TRUE, includes the full fit object in the output (default: FALSE).
 #'
 #' @return A list of data frames summarizing posterior stage probabilities for each model combination.
 #'
 #' @details
 #' Computes posterior means and credible intervals for stage probabilities using HDInterval and dplyr. Results are formatted for each cluster/mixture combination.
 #'
 #' @examples
 #' \dontrun{
 #' # Summarize posterior for a fitted model
 #' summ <- posterior_summary(fit)
 #' }
 #'
 #' @importFrom dplyr bind_cols select filter
 #' @importFrom HDInterval hdi
 #' @export
posterior_summary = function(fit, M = NULL, G = NULL, cred_mass = 0.95, keep_fit = FALSE) {

  if(is.null(M)) M = fit$args$M
  if(is.null(G)) G = fit$args$G

  cluster_dim_list = expand.grid(G = G, M = M)
  cluster_dim_list = lapply(1:nrow(cluster_dim_list), function(i){
    c(G = cluster_dim_list[i, "G"], M = cluster_dim_list[i, "M"])
  })

  model_names = lapply(cluster_dim_list, function(cluster_dim){
    paste0("G=", cluster_dim["G"], ", M=", cluster_dim["M"])
  })

  post_summ_list = lapply(model_names, function(model_name){

    model = fit$models[[model_name]]
    sample = model$sample_list

    M = model$model_info$M %>% as.numeric()
    n_time = fit$model_data$n_time


    prob_mean = sample$stage_prob |> comp_post_stat(mean)
    prob_lower = sample$stage_prob |>
      comp_post_stat(stat_function = function(x){HDInterval::hdi(x)[1]})
    prob_upper = sample$stage_prob |>
      comp_post_stat(stat_function = function(x){HDInterval::hdi(x)[2]})

    time = fit$args$time |> unique() |> rep(times = M)
    w = (1:M) |> rep(each = n_time)

    w_class = model$w_class

    prob_stage_summ = dplyr::bind_cols(
      prob_mean |> trans_summ_prob(col_name = "mean", time = time, w = w),
      prob_lower |> trans_summ_prob(col_name = "li", time = time, w = w) |> dplyr::select(li),
      prob_upper |> trans_summ_prob(col_name = "ui", time = time, w = w) |> dplyr::select(ui)
    ) |>
      dplyr::filter(w %in% unique(w_class))

    z_class = model$z_class

    post_mean = list(
      alpha = sample$alpha |> comp_post_stat(mean),
      mu = sample$mu |> comp_post_stat(mean),
      sigma = sample$sigma |> comp_post_stat(mean)
    )

    out = list(
      w_class = w_class,
      z_class = z_class,
      prob_stage_summ = prob_stage_summ,
      post_mean = post_mean
    )

    out

  })

  names(post_summ_list) = model_names %>% unlist()
  return(post_summ_list)
}

