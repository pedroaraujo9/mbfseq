test_that("function works", {

  M = 3
  z = sim_data$data$true_z
  id = sim_data$data$id
  time = sim_data$data$time
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
  w = NULL

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
    lambda = lambda,
    init_list = init_list,
    epsilon = 1
  )

  alpha = sample_list$alpha[1,,]
  pw = sample_list$pw[1,]
  z = sim_data$data$true_z
  model_data = model_data

  new_w = update_w(
    alpha = alpha,
    pw = pw,
    z = z,
    model_data = model_data
  ) %>%
    expect_no_error() %>%
    expect_no_warning() %>%
    expect_no_message() %>%
    expect_type("list")

})
