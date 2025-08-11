test_that("function works", {

  data = sim_data$data
  time = data$time
  id = data$id
  x = NULL
  z = data$true_z
  G = NULL
  w = NULL
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

  set.seed(1)

  options = list(
    init_points = 5,
    iters_n = 5,
    acq = "ucb",
    kappa = 2.57,
    length_out = 10
  )


  fit1 = calibrate_lambda(
    bounds = c(0.01, 1),
    model_data = model_data,
    method = "bayes",
    fixed_sd = 10,
    options = options
  ) |>
    expect_no_error()

  fit2 = calibrate_lambda(
    bounds = c(0.01, 1),
    model_data = model_data,
    method = "grid",
    fixed_sd = 10,
    options = options
  ) |>
    expect_no_error()

  expect_equal(fit1$best_lambda, fit2$best_lambda, tolerance = 0.1)

})
