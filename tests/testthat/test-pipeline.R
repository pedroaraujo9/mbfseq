test_that("clustering w", {
  cluster_dim = c("M" = 3, "G" = 3)
  z = NULL
  w = sim_data$true_w
  x = sim_data$x
  id = sim_data$data$id
  time = sim_data$data$time
  iters = 80
  burn_in = iters/2
  thin = 2
  lambda = NULL
  n_basis = 10
  init_list = NULL
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
  seed = 1

  fit = pipeline(
    cluster_dim = cluster_dim,
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
    config = config,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  mclust::adjustedRandIndex(fit$fit$w_class, sim_data$true_w) |>
    expect_equal(1)
})

test_that("clustering z", {
  cluster_dim = c("M" = 3, "G" = 3)
  z = NULL
  w = NULL
  x = sim_data$x
  id = sim_data$data$id
  time = sim_data$data$time
  iters = 80
  burn_in = iters/2
  thin = 2
  lambda = NULL
  n_basis = 10
  init_list = NULL
  config = list(
    bounds = c(0.01, 10),
    lambda_start = 1,
    n_points = 20,
    n_start = 10,
    n_start_iters = 10,
    n_start_cores = 1,
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )
  verbose = FALSE
  seed = 1

  fit = pipeline(
    cluster_dim = cluster_dim,
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
    config = config,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  mclust::adjustedRandIndex(fit$fit$w_class, sim_data$true_w) |>
    expect_equal(1)

  mclust::adjustedRandIndex(fit$fit$z_class, sim_data$data$true_z) |>
    expect_gt(0.4)
})







