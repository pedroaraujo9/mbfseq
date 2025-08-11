test_that("function works", {
  x = sim_data$x_random
  z = NULL
  id = sim_data$data$id
  time = sim_data$data$time
  w = NULL
  G = 3
  M = 2
  n_basis = 10
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

  mu = sim_data$true_mu
  sigma = rep(1, G)

  set.seed(1)
  z_up = update_z(mu = mu, sigma = sigma, model_data = model_data) |>
    expect_no_error() |>
    expect_no_message()

  z_up$z
  mclust::adjustedRandIndex(z_up$z, sim_data$data$z_random) |>
    expect_gt(0.5)
})
