test_that("function works", {

  x = sim_data$x_random
  z = sim_data$data$z_random
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

  sigma = rep(1, G)

  mu = update_mu(
    z = z,
    sigma = sigma,
    model_data = model_data,
    mu_fixed_sd = mu_fixed_sd
  ) |>
    expect_no_error() |>
    expect_no_message()

  mean(abs((mu - sim_data$true_mu)/sim_data$true_mu)) %>%
    expect_lt(0.05)
})
