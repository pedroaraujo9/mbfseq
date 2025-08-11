test_that("function works", {

  M = 3
  z = sim_data$data$true_z
  id = sim_data$data$id
  time = sim_data$data$time
  w = sim_data$true_w
  iters = 1000
  thin = 1
  burn_in = iters/2
  lambda = 1
  n_basis = 10
  fixed_sd = 10
  epsilon = 1
  intercept = FALSE
  init_list = NULL
  verbose = FALSE

  model_data = create_model_data(
    id = id,
    time = time,
    z = z,
    w = w,
    n_basis = n_basis,
    G = max(z),
    M = M,
    intercept = intercept
  )

  sample_list = create_sample_list(
    iters = iters,
    model_data = model_data,
    init_list = init_list,
    epsilon = 1
  )

  alpha = sample_list$alpha[1,,]
  lambda = sample_list$lambda[1,,]
  z = sim_data$data$true_z
  w = sim_data$true_w
  fixed_sd = 10
  intercept = FALSE

  update_alpha(
    alpha = alpha,
    lambda = lambda,
    z = z,
    w = w,
    model_data = model_data,
    fixed_sd = fixed_sd
  ) %>%
    expect_no_error() %>%
    expect_no_message()

})



