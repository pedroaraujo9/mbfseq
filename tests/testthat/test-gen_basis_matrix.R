test_that("function works", {

  n_basis = 12
  time = sim_data$data$time
  order = 1

  gen_basis_matrix(
    n_basis = n_basis,
    time = time
  ) |>
    expect_no_error() |>
    expect_no_message()

  #### test 2
  n_basis = 3

  gen_basis_matrix(
    n_basis = n_basis,
    time = time
  ) |>
    expect_no_error() |>
    expect_no_message()

})
