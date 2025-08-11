test_that("funciton works", {
  x = sim_data$x_random
  z = NULL
  id = sim_data$data$id
  time = sim_data$data$time
  w = NULL
  G = 3
  M = 2
  n_basis = 10
  intercept = FALSE
  mu_fixed_sd = 1000

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

  z = sim_data$data$z_random
  mu = sim_data$true_mu
  sigma_a = 1
  sigma_b = 1

  update_sigma(
    mu = mu,
    z = z,
    sigma_a = sigma_a,
    sigma_b = sigma_b,
    model_data = model_data
  ) |>
    expect_no_error() |>
    expect_no_message()

})
