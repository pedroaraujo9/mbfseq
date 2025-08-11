test_that("to estimate model parameters", {

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

  model_data = expect_no_error({
    create_model_data(
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
  })

  expect_equal(model_data$B_expand[, 1:15], model_data$B_expand[, 16:30],
               tolerance = 1e-6
  )

  expect_equal(model_data$B_expand[, 1:15], model_data$B_expand[, 31:45],
               tolerance = 1e-6
  )

  expect_equal(model_data$G, 3)
  expect_equal(model_data$M, 3)
  expect_equal(model_data$n_basis, 15)
  expect_equal(model_data$n_id, 30)
  expect_equal(model_data$n_time, 50)
  expect_equal(dim(model_data$B_unique), c(length(unique(data$time)), 15))
  expect_null(model_data$x)

})

test_that("to estimate w", {

  data = sim_data$data

  time = data$time
  id = data$id
  x = NULL
  z = data$true_z
  w = NULL
  G = NULL
  M = 7
  n_basis = 15
  intercept = FALSE

  model_data = expect_no_error({
    create_model_data(
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
  })

  expect_equal(model_data$G, 3)
  expect_equal(model_data$M, 7)
  expect_null(model_data$w)

})

