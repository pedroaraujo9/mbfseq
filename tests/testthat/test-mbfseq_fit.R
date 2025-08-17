test_that("clustering w", {
  G = 3
  M = c(2, 3)
  z = sim_data$data$true_z
  w = NULL
  x = NULL
  id = sim_data$data$id
  time = sim_data$data$time
  iters = 80
  burn_in = iters/2
  thin = 2
  lambda = NULL
  n_basis = 10
  init_list = NULL
  n_cores = 1
  config = list(
    bounds = c(0.01, 10),
    lambda_start = 1,
    n_points = 20,
    n_start = 5,
    n_start_iters = 10,
    n_start_cores = 1,
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )
  verbose = FALSE
  seed = NULL

  fit = fit_mbfseq(
    G = G,
    M = M,
    z = z,
    w = w,
    x = x,
    id = id,
    time = time,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    lambda = lambda,
    n_basis = n_basis,
    init_list = init_list,
    n_cores = n_cores,
    config = config,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  mclust::adjustedRandIndex(fit$models$`G=3, M=3`$w_class, sim_data$true_w) |>
    expect_equal(1)

})

test_that("clustering z given w", {
  G = 3
  z = NULL
  M = NULL
  w = sim_data$true_w
  x = sim_data$x
  id = sim_data$data$id
  time = sim_data$data$time
  iters = 100
  burn_in = iters/2
  thin = 5
  lambda = NULL
  n_basis = 10
  init_list = NULL
  n_cores = 1
  config = list(
    bounds = c(0.01, 10),
    lambda_start = 1,
    n_points = 20,
    n_start = 5,
    n_start_iters = 10,
    n_start_cores = 1,
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )
  verbose = FALSE
  seed = NULL

  fit = fit_mbfseq(
    G = G,
    M = M,
    z = z,
    w = w,
    x = x,
    id = id,
    time = time,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    lambda = lambda,
    n_basis = n_basis,
    init_list = init_list,
    n_cores = n_cores,
    config = config,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  mclust::adjustedRandIndex(fit$models$`G=3, M=3`$w_class, sim_data$true_w) |>
    expect_equal(1)
})
