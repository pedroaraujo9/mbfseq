test_that("single z clustering works", {
  x = sim_data$x_random
  z = NULL
  id = sim_data$data$id
  time = sim_data$data$time
  w = NULL
  G = 3
  M = 2
  n_basis = 10
  intercept = FALSE
  mu_prior_sd = 1000

  model_data = create_model_data(
    time = time,
    id = id,
    x = x,
    z = z,
    w = w,
    G = G,
    M = M,
    n_basis = n_basis,
    intercept = intercept
  )

  iters = 100
  burn_in = iters/2
  thin = 2
  epsilon = 1
  mu_fixed_sd = 100
  seed = 1
  lambda = 1
  init_list = NULL
  verbose = FALSE
  sigma_a = 1
  sigma_b = 1

  fit = clustz_single_run(
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    model_data = model_data,
    epsilon = epsilon,
    mu_fixed_sd = mu_fixed_sd,
    lambda = lambda,
    init_list = init_list,
    sigma_a = sigma_a,
    sigma_b = sigma_b,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  mclust::adjustedRandIndex(fit$z_class, sim_data$data$z_random) |>
    expect_gt(0.5)
})
