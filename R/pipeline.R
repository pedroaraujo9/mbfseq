pipeline = function(cluster_dim,
                    z,
                    w,
                    x,
                    id,
                    time,
                    iters,
                    burn_in,
                    thin,
                    lambda,
                    n_basis,
                    init_list,
                    config = list(
                      bounds = c(0.01, 10),
                      lambda_start = 1,
                      n_points = 20,
                      n_start = 30,
                      n_start_iters = 20,
                      n_start_cores = 1,
                      epsilon_w = 1,
                      beta_sd = sqrt(10),
                      mu_sd = sqrt(10),
                      sigma_a = 1,
                      sigma_b = 1
                    ),
                    verbose = TRUE,
                    seed = NULL) {

  if(is.null(seed)) seed = sample(1:10000, size = 1)
  set.seed(seed)

  g = cluster_dim["G"]
  m = cluster_dim["M"]

  model_data = create_model_data(
    time = time,
    id = id,
    x = x,
    z = z,
    w = w,
    G = g,
    M = m,
    n_basis = n_basis,
    intercept = FALSE
  )

  # find optimal lambda if lambda is not provided
  if(is.null(lambda)) {

    if(verbose) cat(paste0("G = ", g, ", M = ", m, " - Finding lambda\n"))

    opt_lambda = calibrate_lambda(
      z = z,
      w = w,
      model_data = model_data,
      config = config
    )

  }else{

    opt_lambda = list(best_lambda = lambda)

  }

  # find inits
  if(verbose) cat(paste0("G = ", g, ", M = ", m, " - Finding inits\n"))

  init = find_init(
    n_start = config$n_start,
    iters = config$n_start_iters,
    n_cores = config$n_start_cores,
    model_data = model_data,
    lambda = opt_lambda$best_lambda,
    init_list = init_list,
    priors = config,
    seed = NULL
  )

  # final run
  run = single_run(
    model_data = model_data,
    lambda = opt_lambda$best_lambda,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    init_list = init$init_list,
    priors = config,
    verbose = verbose,
    seed = NULL
  )

  out = list(
    opt_lambda = opt_lambda,
    opt_init = init,
    fit = run,
    seed = seed
  )

  return(out)

}
