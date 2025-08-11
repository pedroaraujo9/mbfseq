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

  iters = 1000
  lambda = 1
  init_list = NULL
  seed = NULL

  sample_list = create_sample_list(
    iters = iters,
    model_data = model_data,
    init_list = init_list
  )

  z = sample_list$z[1, ]
  w = sample_list$w[1, ]
  alpha = sample_list$alpha[1,,]
  mu = sample_list$mu[1,,]
  sigma = sample_list$sigma[1,]

  priors = list(
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )

  eval_seq_logpost(
    z = z,
    alpha = alpha,
    w = w,
    lambda = lambda,
    fixed_sd = priors$beta_sd,
    model_data = model_data
  ) |>
    expect_no_error() |>
    expect_no_message()


  eval_logpost(
    z = z,
    w = w,
    alpha = alpha,
    mu = mu,
    sigma = sigma,
    lambda = lambda,
    priors = priors,
    model_data = model_data
  ) |>
    expect_no_error() |>
    expect_no_message()

})

test_that("clustering only w", {
  data = sim_data$data
  time = data$time
  id = data$id
  x = sim_data$x_random
  z = sim_data$data$true_z
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

  iters = 1000
  lambda = 1
  init_list = NULL
  seed = NULL

  sample_list = create_sample_list(
    iters = iters,
    model_data = model_data,
    init_list = init_list
  )

  z = sample_list$z[1, ]
  w = sample_list$w[1, ]
  alpha = sample_list$alpha[1,,]
  mu = sample_list$mu[1,,]
  sigma = sample_list$sigma[1,]

  priors = list(
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )

  eval_seq_logpost(
    z = z,
    alpha = alpha,
    w = w,
    lambda = lambda,
    fixed_sd = priors$beta_sd,
    model_data = model_data
  ) |>
    expect_no_error() |>
    expect_no_message()


  logpost = eval_logpost(
    z = z,
    w = w,
    alpha = alpha,
    mu = mu,
    sigma = sigma,
    lambda = lambda,
    priors = priors,
    model_data = model_data
  ) |>
    expect_no_error() |>
    expect_no_message()

  unique(logpost) |> length() |> expect_equal(1)
})
