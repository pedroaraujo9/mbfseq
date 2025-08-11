test_that("function works", {

  data = sim_data$data
  time = data$time
  id = data$id
  x = sim_data$x
  w = NULL
  z = NULL
  G = 3
  M = 3
  n_basis = 5
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

  iters = 200
  burn_in = 100
  thin = 2
  init_list = NULL
  priors = list(
    epsilon_w = 1,
    beta_sd = sqrt(10),
    mu_sd = sqrt(10),
    sigma_a = 1,
    sigma_b = 1
  )

  #### single core ####
  verbose = FALSE
  seed = 2
  lambda = 1
  n_start = 5
  iters = 20
  n_cores = 1

  inits = find_init(
    n_start = n_start,
    iters = iters,
    n_cores = n_cores,
    model_data = model_data,
    lambda = lambda,
    thin = thin,
    init_list = init_list,
    priors = priors,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  inits$logpost %>%
    length() %>%
    expect_equal(n_start)

  #### multi core ####
  verbose = FALSE
  seed = 2
  lambda = 1
  n_start = 5
  iters = 20
  n_cores = 2

  inits = find_init(
    n_start = n_start,
    iters = iters,
    n_cores = n_cores,
    model_data = model_data,
    lambda = lambda,
    thin = thin,
    init_list = init_list,
    priors = priors,
    seed = seed
  )

})
