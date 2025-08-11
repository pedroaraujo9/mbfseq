test_that("estimating alpha", {

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

  iters = 1000
  init_list = NULL
  seed = NULL

  sample_list = expect_no_error({
    create_sample_list(
      iters = iters,
      model_data = model_data,
      init_list = init_list
    )
  })

  expect_equal(dim(sample_list$alpha), c(1000, n_basis*M, model_data$G))
  expect_true(
    all(apply(sample_list$w, MARGIN = 1,  function(x) x == w))
  )
  expect_true(
    all(apply(sample_list$z, MARGIN = 1,  function(x) x == z))
  )

})

test_that("estimating w", {
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
  init_list = NULL
  seed = NULL

  sample_list = expect_no_error({
    create_sample_list(
      iters = iters,
      model_data = model_data,
      init_list = init_list
    )
  })

  is.na(sample_list$w[-1, ]) %>%
    all() %>%
    expect_true()

  expect_true(
    all(apply(sample_list$z, MARGIN = 1,  function(x) x == z))
  )

  sample_list$mu %>% expect_null()
  sample_list$sigma %>% expect_null()

})

test_that("estimating w and z", {
  data = sim_data$data

  time = data$time
  id = data$id
  x = sim_data$x
  z = NULL
  w = NULL
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

  iters = 1000
  init_list = NULL
  seed = NULL

  sample_list = expect_no_error({
    create_sample_list(
      iters = iters,
      model_data = model_data,
      init_list = init_list
    )
  })

  is.na(sample_list$w[-1, ]) %>%
    all() %>%
    expect_true()

  is.na(sample_list$z[-1, ]) %>%
    all() %>%
    expect_true()

  expect_false(is.null(sample_list$mu))
  expect_false(is.null(sample_list$sigma))

})

test_that("init values", {
  data = sim_data$data

  time = data$time
  id = data$id
  x = sim_data$x
  z = NULL
  w = NULL
  G = 3
  M = 3
  n_basis = 15
  intercept = FALSE

  w_init = sim_data$true_w
  z_init = data$true_z
  alpha_init = matrix(rnorm(n_basis * M * G), nrow = n_basis * M, ncol = G)
  alpha_init[,G] = 0

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

  init_list = list(
    z = z_init,
    w = w_init,
    alpha = alpha_init
  )
  seed = NULL

  sample_list = expect_no_error({
    create_sample_list(
      iters = iters,
      model_data = model_data,
      init_list = init_list
    )
  })

  expect_equal(sample_list$w[1,], w_init)
  expect_equal(sample_list$z[1,], z_init)
  expect_equal(sample_list$alpha[1,,], alpha_init)

})





