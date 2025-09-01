#' Fit Sequential Model
#'
#' Fits a multi-basis function sequential model for given combinations of
#' `G` and `M`, performing lambda calibration, initialization search,
#' and final model fitting. Supports parallel processing for model combinations.
#'
#' @param G Numeric vector of possible values for `G`.
#' @param M Numeric vector of possible values for `M`.
#' @param z Optional sequence values to be clustered.
#' @param w Optional sequence clusters.
#' @param x Optional covariate matrix to find `z`.
#' @param id Vector of subject or group identifiers.
#' @param time Vector of time points corresponding to the observations.
#' @param iters Integer. Total number of MCMC iterations. Default is 1000.
#' @param burn_in Integer. Number of burn-in iterations. Default is `iters/2`.
#' @param thin Integer. Thinning interval for MCMC sampling. Default is 5.
#' @param lambda Optional numeric value for the regularization parameter. If `NULL`,
#'   the optimal lambda is estimated via [calibrate_lambda()].
#' @param n_basis Integer. Number of spline basis functions (includes intercept). Default is 10.
#' @param init_list Optional list of initial values for model parameters.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing. Default is 1.
#' @param config List of additional configuration parameters:
#'   - `single_group`: Logical. If `TRUE`, assumes a single group for clustering. Default is `FALSE`.
#'   - `n_start`: Number of random initializations.
#'   - `n_start_iter`: Number of iterations for initialization search.
#'   - `n_start_cores`: Number of cores for initialization search.
#'   - `epsilon_w`, `beta_sd`, `mu_sd`, `sigma_a`, `sigma_b`: Priors and scaling constants.
#'   - `bounds`: Range for lambda search.
#'   - `n_grid`: Number of grid points for lambda search.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default is TRUE.
#' @param seed Optional integer seed for reproducibility. If `NULL`, a random seed is generated.
#'
#' @return A list containing:
#'   - `fit`: List of fitted model objects.
#'   - `init`: List of initialization results.
#'   - `lambda_opt`: List of optimal lambda search results.
#'   - `best_lambda`: Numeric vector of best lambda values.
#'   - `run_time`: Total runtime as a `difftime` object.
#'   - `seed`: Random seed used.
#'   - `metrics`: Formatted model performance metrics.
#'   - `model_data`: Minimal model data object.
#'
#' @details
#' This function iterates over all combinations of `G` and `M` values,
#' fits the corresponding model using [single_run()], and computes performance metrics
#' with [comp_metrics()].
#'
#' If `lambda` is `NULL`, [calibrate_lambda()] is used to find the optimal value.
#' Initialization values are obtained using [find_init()].
#'
#' @examples
#' \dontrun{
#' result <- fit_mbfseq(G = c(2, 3), M = c(4, 5),
#'                      id = rep(1:10, each = 5),
#'                      time = rep(1:5, times = 10),
#'                      iters = 500, n_cores = 2)
#' }
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom purrr map
#'
#' @export
fit_mbfseq = function(G = NULL,
                      M = NULL,
                      z = NULL,
                      w = NULL,
                      x = NULL,
                      id,
                      time,
                      iters = 1000,
                      burn_in = iters/2,
                      thin = 5,
                      lambda = NULL,
                      n_basis = 10,
                      init_list = NULL,
                      n_cores = 1,
                      config = list(
                        single_group = FALSE,
                        bounds = c(0.01, 10),
                        lambda_start = 1,
                        n_points = 20,
                        n_start = 30,
                        n_start_iters = 20,
                        n_start_cores = 1,
                        epsilon_w = 1,
                        beta_sd = sqrt(10),
                        mu_sd = sqrt(10),
                        sigma_a = 1,
                        sigma_b = 1
                      ),
                      verbose = TRUE,
                      seed = NULL
){

  # Save all argument values in a list
  args = list(
    G = G,
    M = M,
    z = z,
    w = w,
    x = x,
    id = id,
    time = time,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    lambda = lambda,
    n_basis = n_basis,
    init_list = init_list,
    n_cores = n_cores,
    config = config,
    verbose = verbose,
    seed = seed
  )

  if(is.null(seed)) seed = sample(1:10000, size = 1)
  set.seed(seed)

  if(!is.null(z)) {
    G = length(unique(z))
  }

  if(!is.null(w)) {
    M = length(unique(w))
  }

  if(is.factor(z)) {

    z_levels = z %>% levels()
    z = z %>% as.integer()

  }else{

    z_levels = unique(z) %>% sort()

  }

  model_data_min = create_model_data(
    time = time,
    id = id,
    x = x,
    z = z,
    w = w,
    G = min(G),
    M = min(M),
    n_basis = n_basis,
    intercept = FALSE
  ) %>%
    create_model_data_min(M = M, G = G)

  model_data_min$z_levels = z_levels


  # config parallel session
  if(n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  }else{
    future::plan(future::sequential, split = TRUE)
  }

  # list with all models combination
  cluster_dim_list = expand.grid(G = G, M = M)
  cluster_dim_list = lapply(1:nrow(cluster_dim_list), function(i){
    c(G = cluster_dim_list[i, "G"], M = cluster_dim_list[i, "M"])
  })

  # fun model
  init_time = Sys.time()
  runs = future.apply::future_lapply(cluster_dim_list, function(cluster_dim){

  pipeline(
    cluster_dim = cluster_dim,
    z = z,
    w = w,
    x = x,
    id = id,
    time = time,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    lambda = lambda,
    n_basis = n_basis,
    init_list = init_list,
    config = config,
    verbose = verbose,
    seed = NULL
  )

  }, future.seed = TRUE)
  end_time = Sys.time()
  run_time = end_time - init_time

  names(runs) = lapply(cluster_dim_list, function(cluster_dim){
    paste0("G=", cluster_dim["G"], ", M=", cluster_dim["M"])
  }) %>% do.call(rbind, .)

  models = purrr::map(runs, ~{.x$fit})
  lambda_opt = purrr::map(runs, ~{.x$opt_lambda})
  best_lambda = purrr::map(lambda_opt, ~{.x$best_lambda}) %>% format_best_lambda()
  init_opt = purrr::map(runs, ~{.x$opt_init})
  metrics = comp_metrics(models, model_data_min = model_data_min) %>% format_metrics()

  out = list(
    models = models,
    init = init_opt,
    lambda_opt = lambda_opt,
    best_lambda = best_lambda,
    run_time = run_time,
    seed = seed,
    metrics = metrics,
    model_data = model_data_min,
    args = args
  )

  return(out)

}
