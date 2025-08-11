test_that("data", {

  df = sim_data$data

  Z = expect_no_error({
    gen_dummy(df$true_z, n_cat = 3, intercept = FALSE)
  })

  Zi = expect_no_error({
    gen_dummy(df$true_z, n_cat = 3, intercept = TRUE)
  })

  expect_false(any(is.na(Z)))
  expect_false(any(is.na(Zi)))

})

test_that("indices", {
  expect_equal(
    gen_notpen_index(n_basis = 10, M = 2, order = 2), c(1, 2, 11, 12)
  )

  expect_equal(
    gen_notpen_index(n_basis = 5, M = 1, order = 2), c(1, 2)
  )

  expect_equal(
    gen_notpen_index(n_basis = 5, M = 3, order = 1), c(1, 6, 11)
  )

  expect_equal(
    gen_basis_index(n_basis = 3, M = 2), list(1:3, 4:6)
  )

  expect_equal(
    gen_basis_index(n_basis = 10, M = 3), list(1:10, 11:20, 21:30)
  )

  expect_equal(
    gen_basis_index(n_basis = 12, M = 2), list(1:12, 13:24)
  )

})


test_that("generating precision matrix", {
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

  iters = 1000
  epsilon = 1
  lambda = 2
  init_list = NULL
  seed = 1
  fixed_sd = 10

  sample_list = create_sample_list(
    model_data = model_data,
    iters = iters,
    epsilon = epsilon,
    lambda = lambda,
    init_list = init_list,
    seed = seed
  )

  lambda = sample_list$lambda[1,,]

  inv_cov = gen_inv_cov(
    lambda = lambda,
    model_data = model_data,
    fixed_sd = fixed_sd
  ) |>
    expect_no_error() |>
    expect_no_message()

  expect_equal(length(inv_cov), model_data$G - 1)
  expect_equal(inv_cov[[1]], inv_cov[[2]])

  diag_inv_cov = inv_cov[[1]] |> diag()

  expect_equal(
    diag_inv_cov[1:(n_basis)],
    diag_inv_cov[(n_basis + 1):(n_basis + n_basis)]
  )

  expect_equal(1/diag_inv_cov[model_data$notpen_index], rep(fixed_sd^2, M))


  time = data$time
  id = data$id
  z = sample(1:2, size = length(time), replace = T)

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
  epsilon = 1
  lambda = 2
  init_list = NULL
  seed = 1
  fixed_sd = 10

  sample_list = create_sample_list(
    model_data = model_data,
    iters = iters,
    epsilon = epsilon,
    lambda = lambda,
    init_list = init_list,
    seed = seed
  )

  lambda = sample_list$lambda[1,,]

  inv_cov = gen_inv_cov(
    lambda = lambda,
    model_data = model_data,
    fixed_sd = fixed_sd
  ) |>
    expect_no_error() |>
    expect_no_message()

})


test_that("reg matrix", {

  data = sim_data$data

  time = data$time
  id = data$id
  x = NULL
  w = NULL
  z = data$true_z
  G = NULL
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

  w_vec = data$true_w
  B_expand = model_data$B_expand

  X3 = expect_no_error({
    gen_reg_matrix(
      w_vec = w_vec,
      B_expand = model_data$B_expand,
      M = M
    )
  })

})
