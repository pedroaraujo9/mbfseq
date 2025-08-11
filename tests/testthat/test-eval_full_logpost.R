test_that("function works", {

  data = sim_data$data

  time = data$time
  id = data$id
  x = sim_data$x_random
  z = NULL
  w = sim_data$true_w
  G = 3
  M = 3
  n_basis = 15
  intercept = FALSE

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

  model_data = model_data
  iters = 1000
  epsilon = 1
  lambda = 10
  init_list = NULL
  seed = NULL

  sample_list = create_sample_list(
    iters = iters,
    epsilon = epsilon,
    model_data = model_data,
    lambda = lambda,
    init_list = init_list
  )

  z = sample_list$z[1, ]
  w = sample_list$w[1, ]
  alpha = sample_list$alpha[1,,]
  mu = sample_list$mu[1,,]
  sigma = sample_list$sigma[1,]
  lambda = sample_list$lambda[1,,]
  sigma_a = sigma_b = 1
  spline_fixed_sd = 10
  mu_fixed_sd = 10

  eval_full_logpost(
    z = z,
    w = w,
    alpha = alpha,
    mu = mu,
    sigma = sigma,
    lambda = lambda,
    spline_fixed_sd = spline_fixed_sd,
    mu_fixed_sd = mu_fixed_sd,
    sigma_a = sigma_a,
    sigma_b = sigma_b,
    model_data = model_data
  ) |>
    expect_no_error() |>
    expect_no_message()

})
