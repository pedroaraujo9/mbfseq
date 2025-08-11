test_that("function works", {

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

  n_start = 10
  max_iter_start = 10
  n_cores_init = 1
  set.seed(1)

  find_init = clustseq_find_init(
    model_data = model_data,
    n_start = n_start,
    max_iter_start = max_iter_start,
    n_cores_init = n_cores_init,
    w = w,
    lambda = lambda,
    n_basis = n_basis,
    fixed_sd = fixed_sd,
    epsilon = epsilon,
    intercept = intercept,
    init_list = init_list ,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

})

test_that("parallel run", {

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

  n_start = 10
  max_iter_start = 10
  n_cores_init = 5

  find_init = clustseq_find_init(
    model_data = model_data,
    n_start = n_start,
    max_iter_start = max_iter_start,
    n_cores_init = n_cores_init,
    w = w,
    lambda = lambda,
    n_basis = n_basis,
    fixed_sd = fixed_sd,
    epsilon = epsilon,
    intercept = intercept,
    init_list = init_list ,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

})

