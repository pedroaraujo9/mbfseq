test_that("find_map works", {

  time = sim_data$data$time
  id = sim_data$data$id
  x = NULL
  z = sim_data$data$true_z
  w = sim_data$true_w
  G = NULL
  M = NULL
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

  model_path = system.file(
    "extdata", "model-stan.rds", package = "mbfseq"
  )

  stan_model = readRDS(model_path)

  map = find_map(
    z = z,
    w = w,
    lambda = 2,
    stan_model = stan_model,
    model_data = model_data
  ) %>%
    expect_no_error() %>%
    expect_no_message()

})

test_that("z and w known", {

  config = list(
    bounds = c(0.01, 10),
    n_points = 10,
    n_start_iters = NULL,
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )

  z = sim_data$data$true_z
  w = sim_data$true_w
  time = sim_data$data$time
  id = sim_data$data$id
  x = NULL
  G = NULL
  M = NULL
  n_basis = 10
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

  lambda_opt1 = calibrate_lambda(
    z = z,
    w = w,
    model_data = model_data,
    config = config
  ) |> expect_no_error() |>
    expect_no_message()

})

test_that("z known and w unknown", {

  config = list(
    bounds = c(0.01, 10),
    n_points = 10,
    n_start_iters = NULL,
    start_lambda = NULL,
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )

  z = sim_data$data$true_z
  w = NULL
  time = sim_data$data$time
  id = sim_data$data$id
  x = NULL
  G = NULL
  M = 3
  n_basis = 10
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

  lambda_opt2 = calibrate_lambda(
    z = z,
    w = w,
    model_data = model_data,
    config = config
  ) |> expect_no_error() |>
    expect_no_message()

})

test_that("z unknown and w known", {

  config = list(
    bounds = c(0.01, 5),
    n_points = 10,
    n_start_iter = 10,
    lambda_start = 1,
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )

  z = NULL
  w = sim_data$true_w
  time = sim_data$data$time
  id = sim_data$data$id
  x = sim_data$x
  G = 3
  M = NULL
  n_basis = 10
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

  lambda_opt3 = calibrate_lambda(
    z = z,
    w = w,
    model_data = model_data,
    config = config
  ) |> expect_no_error() |>
    expect_no_message()

  lambda_opt3$best_lambda

})


test_that("both z unknown and w known", {

  config = list(
    bounds = c(0.01, 5),
    n_points = 10,
    n_start_iters = 10,
    lambda_start = 1,
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )

  z = NULL
  w = NULL
  time = sim_data$data$time
  id = sim_data$data$id
  x = sim_data$x
  G = 3
  M = 3
  n_basis = 10
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

  lambda_opt3 = calibrate_lambda(
    z = z,
    w = w,
    model_data = model_data,
    config = config
  ) |> expect_no_error() |>
    expect_no_message()

})
