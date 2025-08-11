test_that("function works", {
  set.seed(1)

  M = 3
  z = sim_data$data$true_z
  id = sim_data$data$id
  time = sim_data$data$time
  w = NULL
  iters = 100
  thin = 1
  burn_in = iters/2
  lambda = 1
  n_basis = 10
  fixed_sd = 10
  epsilon = 1
  intercept = FALSE
  init_list = NULL
  verbose = FALSE
  seed = NULL

  max_iter_start = 10
  n_start = 10
  n_cores = 1
  n_cores_init = 1
  bounds = c(0.01, 2)

  fit = cluster_seq(
    M = M,
    z = z,
    id = id,
    time = time,
    w = w,
    lambda = lambda,
    n_basis = n_basis,
    fixed_sd = fixed_sd,
    iters = iters,
    thin = thin,
    burn_in = burn_in,
    epsilon = epsilon,
    max_iter_start = max_iter_start,
    n_start = n_start,
    n_cores = n_cores,
    n_cores_init = n_cores_init,
    intercept = intercept,
    bounds = bounds,
    verbose = FALSE
  ) |>
    expect_no_error()

  class_est = fit$fits$M3$w_class
  expect_equal(mclust::adjustedRandIndex(class_est, sim_data$true_w), 1)
})


test_that("parallell", {
  set.seed(1)

  M = c(2, 3, 4)
  z = sim_data$data$true_z
  id = sim_data$data$id
  time = sim_data$data$time
  w = NULL
  iters = 100
  thin = 1
  burn_in = iters/2
  lambda = 1
  n_basis = 10
  fixed_sd = 10
  epsilon = 1
  intercept = FALSE
  init_list = NULL
  verbose = FALSE
  seed = NULL

  max_iter_start = 10
  n_start = 10
  n_cores = 3
  n_cores_init = 1

  fit = cluster_seq(
    M = M,
    z = z,
    id = id,
    time = time,
    w = w,
    lambda = lambda,
    n_basis = n_basis,
    fixed_sd = fixed_sd,
    iters = iters,
    thin = thin,
    burn_in = burn_in,
    epsilon = epsilon,
    max_iter_start = max_iter_start,
    n_start = n_start,
    n_cores = n_cores,
    n_cores_init = n_cores_init,
    intercept = intercept,
    verbose = FALSE
  ) |>
    expect_no_error()
})
