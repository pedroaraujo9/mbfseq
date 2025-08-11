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

  model_data = create_model_data(
    time = sim_data$data$time,
    id = sim_data$data$id,
    z = sim_data$data$true_z,
    M = M,
    n_basis = 10
  )

  G = model_data$G
  M = model_data$M
  B_expand = model_data$B_expand
  id = model_data$id_time_df$id
  basis_index = model_data$basis_index
  n_basis = model_data$n_basis
  B_unique = model_data$B_unique
  notpen_index = model_data$notpen_index

  sample_list = create_sample_list(
    iters = iters,
    model_data = model_data,
    lambda = lambda,
    init_list = init_list,
    epsilon = 1
  )

  alpha = sample_list$alpha[1,,]
  lambda = sample_list$lambda[1,,]

  inv_cov_list = gen_inv_cov(
    lambda,
    model_data = model_data
  )

  W = gen_dummy(w, M, intercept = intercept)
  X = kronecker(W, B_unique)
  Z = gen_dummy(z, n_cat = G, intercept = FALSE)

  g = 1

  linear_pred = X %*% alpha

  augment_alpha(
    g = g,
    Z = Z,
    X = X,
    linear_pred = linear_pred,
    inv_cov_list = inv_cov_list
  ) |>
    expect_no_error() |>
    expect_no_message()

})


