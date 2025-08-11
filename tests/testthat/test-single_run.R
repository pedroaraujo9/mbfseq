test_that("without clustering", {

  data = sim_data$data
  time = data$time
  id = data$id
  x = NULL
  w = sim_data$true_w
  z = data$true_z
  G = NULL
  M = NULL
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
  verbose = FALSE
  seed = 1
  lambda = 1

  fit = single_run(
    model_data = model_data,
    lambda = lambda,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    init_list = init_list,
    priors = priors,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  est_stage_prob= fit$sample_list$stage_prob %>% comp_post_stat()
  est_stage_prob %>% dim()

  true_stage_prob = sim_data$probs %>%
    spread(cat, prob) %>%
    arrange(w, time) %>%
    select(`1`:`3`) %>%
    as.matrix()

  cor(as.numeric(true_stage_prob), as.numeric(est_stage_prob)) %>%
    expect_gt(0.9)
})

test_that("clustering w", {
  data = sim_data$data
  time = data$time
  id = data$id
  x = NULL
  w = NULL
  z = data$true_z
  G = NULL
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
  verbose = FALSE
  seed = 1
  lambda = 1

  fit = single_run(
    model_data = model_data,
    lambda = lambda,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    init_list = init_list,
    priors = priors,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  mclust::adjustedRandIndex(fit$w_class, sim_data$true_w) |>
    expect_equal(1)
})

test_that("clustering z and w", {
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
  verbose = FALSE
  seed = 1
  lambda = 1

  fit = single_run(
    model_data = model_data,
    lambda = lambda,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    init_list = init_list,
    priors = priors,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  w_class = fit$w_class
  z_class = fit$z_class

  mclust::adjustedRandIndex(w_class, sim_data$true_w) |>
    expect_equal(1)

  mclust::adjustedRandIndex(z_class, data$true_z) |>
    expect_gt(0.4)
})


