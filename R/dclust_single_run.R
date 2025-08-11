dclust_single_run = function(iters = 1000,
                             burn_in = iters/2,
                             thin = 10,
                             model_data,
                             epsilon = 1,
                             mu_fixed_sd = sqrt(10),
                             spline_fixed_sd = sqrt(10),
                             lambda = 1,
                             sigma_a = 1,
                             sigma_b = 1,
                             init_list = NULL,
                             verbose = FALSE,
                             intercept = FALSE,
                             seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  M = model_data$M
  G = model_data$G

  sample_list = create_sample_list(
    iters = iters,
    model_data = model_data,
    lambda = lambda,
    init_list = init_list,
    epsilon = epsilon,
    seed = NULL
  )

  print_h = floor(iters/10)

  fit = kmeans(model_data$x, centers = G)
  mu_init = fit$centers %>% as.matrix()
  z_init = fit$cluster

  sample_list$z[1, ] = z_init
  sample_list$mu[1,,] = mu_init
  sample_list$w[1,] = model_data$z_ham_dist |> hclust(method = "ward.D") |> cutree(k = M)
  names(sample_list$w[1,]) = model_data$id_unique

  logp = matrix(nrow = iters, ncol = 6)
  logp_init = eval_full_logpost(
    z = sample_list$z[1, ],
    w = sample_list$w[1, ],
    alpha = sample_list$alpha[1,,],
    mu = sample_list$mu[1,,],
    sigma = sample_list$sigma[1, ],
    model_data = model_data,
    lambda = sample_list$lambda[1,,],
    sigma_a = sigma_a,
    sigma_b = sigma_b,
    spline_fixed_sd = spline_fixed_sd,
    mu_fixed_sd = mu_fixed_sd
  )

  logp[1, ] = logp_init
  colnames(logp) = names(logp_init)

  for(i in 2:iters) {

    if(verbose) {
      if (i %% print_h == 0) {
        cat(sprintf("M = %d - Iteration %d / %d\n", M, i, iters))
        utils::flush.console()
      }
    }

    #### update alpha, w, pw ####
    sample_list$alpha[i,,] = update_alpha(
      alpha = sample_list$alpha[i-1,,],
      lambda = sample_list$lambda[i-1,,],
      z = sample_list$z[i-1, ],
      w = sample_list$w[i-1, ],
      model_data = model_data,
      fixed_sd = spline_fixed_sd,
      intercept = intercept
    )

    sample_list$prob[i,,] = comp_prob(
      alpha = sample_list$alpha[i,,],
      model_data = model_data
    )

    sample_list$pw[i, ] = update_pw(
      w = sample_list$w[i-1,],
      epsilon = epsilon,
      model_data = model_data
    )

    new_w = update_w(
      alpha = sample_list$alpha[i,,],
      pw = sample_list$pw[i, ],
      z = sample_list$z[i-1, ],
      model_data = model_data
    )

    sample_list$w[i, ] = new_w$w
    sample_list$id_marg_logp[i,] = new_w$marg_log_like
    sample_list$w_post_prob[i,,] = new_w$w_post_prob

    #### update z, mu, sigma ####
    sample_list$mu[i,,] = update_mu(
      z = sample_list$z[i-1, ],
      sigma = sample_list$sigma[i-1, ],
      model_data = model_data,
      mu_fixed_sd = mu_fixed_sd
    )

    sample_list$sigma[i, ] = update_sigma(
      mu = sample_list$mu[i,,],
      z = sample_list$z[i-1, ],
      sigma_a = sigma_a,
      sigma_b = sigma_b,
      model_data = model_data
    )

    z_up = update_z(
      mu = sample_list$mu[i,,],
      sigma = sample_list$sigma[i, ],
      w = sample_list$w[i,],
      alpha = sample_list$alpha[i,,],
      model_data = model_data
    )

    sample_list$z[i, ] = z_up$z
    sample_list$z_post_prob[i,,] = z_up$z_post_prob

    logp[i, ] = eval_full_logpost(
      z = sample_list$z[i, ],
      w = sample_list$w[i, ],
      alpha = sample_list$alpha[i,,],
      mu = sample_list$mu[i,,],
      sigma = sample_list$sigma[i, ],
      model_data = model_data,
      lambda = sample_list$lambda[i,,],
      sigma_a = sigma_a,
      sigma_b = sigma_b,
      spline_fixed_sd = spline_fixed_sd,
      mu_fixed_sd = mu_fixed_sd
    )
  }

  sample_list$logp = logp[seq(burn_in, iters, thin), ]

  if(verbose) {
    cat("\n")
  }

  sample_list = filter_chain(
    sample_list,
    burn_in = burn_in,
    iters = iters,
    thin = thin
  )

  out = list(
    z_class = sample_list$z |> comp_class(),
    w_class = sample_list$w |> comp_class(),
    sample_list = sample_list,
    model_data = model_data
  )

  class(out) = "dclust"
  return(out)
}
