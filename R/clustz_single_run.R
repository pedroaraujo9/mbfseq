clustz_single_run = function(iters = 1000,
                             burn_in = iters/2,
                             thin = 10,
                             model_data,
                             epsilon = 1,
                             mu_fixed_sd = 100,
                             lambda = 1,
                             sigma_a = 1,
                             sigma_b = 1,
                             init_list = NULL,
                             verbose = FALSE,
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

  for(i in 2:iters) {

    if(verbose) {
      if (i %% print_h == 0) {
        cat(sprintf("M = %d - Iteration %d / %d\n", M, i, iters))
        utils::flush.console()
      }
    }

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
      sigma = sample_list$sigma[i,],
      model_data = model_data
    )

    sample_list$z[i, ] = z_up$z
    sample_list$z_post_prob[i,,] = z_up$z_post_prob

  }

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
    sample_list = sample_list,
    model_data = model_data
  )

  class(out) = "dclust"
  return(out)
}
