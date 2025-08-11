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




