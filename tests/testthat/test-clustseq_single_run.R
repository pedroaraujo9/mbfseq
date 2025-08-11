test_that("function works", {

  M = 3
  z = sim_data$data$true_z
  id = sim_data$data$id
  time = sim_data$data$time
  w = sim_data$true_w
  iters = 100
  thin = 1
  burn_in = iters/2
  lambda = 1
  n_basis = 10
  fixed_sd = 10
  epsilon = 1
  intercept = FALSE
  init_list = NULL
  verbose = FALSE
  seed = NULL

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

  set.seed(1)
  fit = clustseq_single_run(
    model_data = model_data,
    w = w,
    iters = iters,
    thin = thin,
    burn_in = burn_in,
    lambda = lambda,
    n_basis = n_basis,
    fixed_sd = fixed_sd,
    epsilon = epsilon,
    intercept = intercept,
    init_list = init_list,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  probs = fit$sample_list$prob |>
    comp_post_stat(mean) |>
    as.data.frame() |>
    dplyr::mutate(time = rep(unique(sim_data$data$time), times = 3),
                  w = rep(1:3, each = 50)) |>
    dplyr::distinct() |>
    tidyr::gather(cat, prob_est, -time, -w) |>
    dplyr::mutate(cat = str_replace(cat, "V", "")) |>
    dplyr::left_join(sim_data$probs, by = c("time", "w", "cat"))

  expect_equal(probs$prob, probs$prob_est, tolerance = 0.15)

  #### posterior summary ####
  post_summ = fit |>
    posterior_summary() |>
    expect_no_error() |>
    expect_no_message()

  post_summ |>
    plot_probability() |>
    expect_no_error() |>
    expect_no_message() |>
    expect_s3_class("ggplot")

})

test_that("function works", {

  M = 3
  z = sim_data$data$true_z
  id = sim_data$data$id
  time = sim_data$data$time
  w = NULL
  iters = 100
  thin = 1
  burn_in = iters/2
  lambda = 1
  n_basis = 10
  fixed_sd = 10
  epsilon = 1
  intercept = FALSE
  init_list = NULL
  verbose = FALSE
  seed = NULL

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

  set.seed(10)
  fit = clustseq_single_run(
    model_data = model_data,
    w = w,
    iters = iters,
    thin = thin,
    burn_in = burn_in,
    lambda = lambda,
    n_basis = n_basis,
    fixed_sd = fixed_sd,
    epsilon = epsilon,
    intercept = intercept,
    init_list = init_list,
    verbose = verbose,
    seed = seed
  ) |>
    expect_no_error() |>
    expect_no_message()

  expect_equal(mclust::adjustedRandIndex(fit$w_class, sim_data$true_w), 1)

})


