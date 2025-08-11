test_that("sample list to estimate w", {

  data = sim_data$data

  time = data$time
  id = data$id
  x = NULL
  z = data$true_z
  w = sim_data$true_w
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

  model_data = model_data
  iters = 1000
  epsilon = 1
  lambda = 10
  init_list = NULL
  seed = NULL

  sample_list = expect_no_error({
    create_sample_list(
      iters = iters,
      epsilon = epsilon,
      model_data = model_data,
      lambda = lambda,
      init_list = init_list
    )
  })

  expect_equal(dim(sample_list$alpha), c(1000, 15*3, 3))
  expect_equal(dim(sample_list$lambda), c(1000, 3, 2))
  expect_equal(dim(sample_list$w), c(1000, model_data$n_id))
  expect_equal(unique(sample_list$lambda[10,3,]), lambda)

})

test_that("inits", {
  data = sim_data$data
  time = data$time
  id = data$id
  x = NULL
  z = data$true_z
  w = NULL
  G = NULL
  M = 3
  n_basis = 15
  intercept = FALSE
  iters = 100
  lambda = 1

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
  alpha_init = matrix(
    rnorm(model_data$n_basis * model_data$M * model_data$G, sd = 0.01),
    nrow = model_data$n_basis * model_data$M, ncol = model_data$G
  )

  w_init = sim_data$true_w
  pw_init = w_init |> table() |> prop.table() |> as.numeric()
  names(pw_init) = NULL
  z_init = NULL

  init_list = list(
    alpha_init = alpha_init,
    w_init = w_init,
    pw_init = pw_init
  )

  sample_list = expect_no_error({
    create_sample_list(
      iters = iters,
      epsilon = epsilon,
      model_data = model_data,
      lambda = lambda,
      init_list = init_list
    )
  })

  expect_equal(sample_list$w[1, ], init_list$w)
  expect_equal(sample_list$pw[1, ], init_list$pw, tolerance = 0.001)
  expect_equal(sample_list$alpha[1, ,1:2], init_list$alpha[, 1:2], tolerance = 0.001)

})


test_that("estimating z", {
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

  model_data = model_data
  iters = 1000
  epsilon = 1
  lambda = 10
  init_list = NULL
  seed = NULL

  sample_list = expect_no_error({
    create_sample_list(
      iters = iters,
      epsilon = epsilon,
      model_data = model_data,
      lambda = lambda,
      init_list = init_list
    )
  }) |>
    expect_no_error() |>
    expect_no_message()

  expect_true(!is.null(sample_list$sigma))

})
