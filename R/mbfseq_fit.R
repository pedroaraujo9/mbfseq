#' Fit Multi-Basis Function Sequential Model
#'
#' Fits a multi-basis function sequential model for given combinations of
#' \code{G} and \code{M}, performing lambda calibration, initialization search,
#' and final model fitting. Supports parallel processing for model combinations.
#'
#' @param G Numeric vector of possible values for \code{G}.
#' @param M Numeric vector of possible values for \code{M}.
#' @param z Optional sequence values to be clustered.
#' @param w Optional sequence clusters.
#' @param x Optional covariate matrix to find z.
#' @param id Vector of subject or group identifiers.
#' @param time Vector of time points corresponding to the observations.
#' @param iters Integer. Total number of MCMC iterations. Default is 1000.
#' @param burn_in Integer. Number of burn-in iterations. Default is \code{iters/2}.
#' @param thin Integer. Thinning interval for MCMC sampling. Default is 5.
#' @param lambda Optional numeric value for the regularization parameter. If \code{NULL},
#'   the optimal lambda is estimated via \code{\link{calibrate_lambda}}.
#' @param n_basis Integer. Number of spline basis functions (includes intercept). Default is 10.
#' @param init_list Optional list of initial values for model parameters.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing. Default is 1.
#' @param config List of additional configuration parameters:
#'   \describe{
#'     \item{\code{n_start}}{Number of random initializations.}
#'     \item{\code{n_start_iter}}{Number of iterations for initialization search.}
#'     \item{\code{n_start_cores}}{Number of cores for initialization search.}
#'     \item{\code{epsilon_w}, \code{beta_sd}, \code{mu_sd}, \code{sigma_a}, \code{sigma_b}}{Priors and scaling constants.}
#'     \item{\code{bounds}}{Range for lambda search.}
#'     \item{\code{n_grid}}{Number of grid points for lambda search.}
#'   }
#' @param verbose Logical. If \code{TRUE}, prints progress messages. Default is TRUE.
#' @param seed Optional integer seed for reproducibility. If \code{NULL}, a random seed is generated.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{fit}}{List of fitted model objects.}
#'   \item{\code{init}}{List of initialization results.}
#'   \item{\code{lambda_opt}}{List of optimal lambda search results.}
#'   \item{\code{best_lambda}}{Numeric vector of best lambda values.}
#'   \item{\code{run_time}}{Total runtime as a \code{difftime} object.}
#'   \item{\code{seed}}{Random seed used.}
#'   \item{\code{metrics}}{Formatted model performance metrics.}
#'   \item{\code{model_data}}{Minimal model data object.}
#' }
#'
#' @details
#' This function iterates over all combinations of \code{G} and \code{M} values,
#' fits the corresponding model using \code{\link{single_run}}, and computes performance metrics
#' with \code{\link{comp_metrics}}.
#'
#' If \code{lambda} is \code{NULL}, \code{\link{calibrate_lambda}} is used to find the optimal value.
#' Initialization values are obtained using \code{\link{find_init}}.
#'
#' @examples
#' \dontrun{
#' result <- mbfseq_fit(G = c(2, 3), M = c(4, 5),
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
mbfseq_fit = function(G = NULL,
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
                        n_start = 30,
                        n_start_iter = 20,
                        n_start_cores = 1,
                        epsilon_w = 1,
                        beta_sd = sqrt(10),
                        mu_sd = sqrt(10),
                        sigma_a = 1,
                        sigma_b = 1,
                        bounds = c(0.01, 5),
                        n_grid = 30
                      ),
                      verbose = TRUE,
                      seed = NULL
){
  if(is.null(seed)) seed = sample(1:10000, size = 1)
  set.seed(seed)

  if(!is.null(z)) {
    G = length(unique(z))
  }

  if(!is.null(w)) {
    M = length(unique(w))
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


  # config parallel session
  if(n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  }else{
    future::plan(future::sequential, split = TRUE)
  }

  # list with all models combination
  models = expand.grid(G = G, M = M)
  models = lapply(1:nrow(models), function(i){
    c(G = models[i, "G"], M = models[i, "M"])
  })

  # fun model
  init_time = Sys.time()
  runs = future.apply::future_lapply(models, function(model){

    g = model["G"]
    m = model["M"]

    model_data = create_model_data(
      time = time,
      id = id,
      x = x,
      z = z,
      w = w,
      G = g,
      M = m,
      n_basis = n_basis,
      intercept = FALSE
    )

    # find optimal lambda
    if(is.null(lambda)) {

      if(verbose) cat(paste0("G = ", g, ", M = ", m, " - Finding lambda\n"))

      opt_lambda = calibrate_lambda(
        bounds = config$bounds,
        model_data = model_data,
        fixed_sd = config$beta_sd,
        method = "grid",
        options = list(length_out = config$n_grid)
      )

      lambda = opt_lambda$best_lambda
    }else{

      opt_lambda = list(best_lambda = lambda)

    }

    # find inits
    if(verbose) cat(paste0("G = ", g, ", M = ", m, " - Finding inits\n"))

    init = find_init(
      n_start = config$n_start,
      iters = config$n_start_iter,
      n_cores = config$n_start_cores,
      model_data = model_data,
      lambda = lambda,
      init_list = init_list,
      priors = config,
      seed = NULL
    )

    # final run
    run = single_run(
      model_data = model_data,
      lambda = lambda,
      iters = iters,
      burn_in = burn_in,
      thin = thin,
      init_list = init$init_list,
      priors = config,
      verbose = verbose,
      seed = NULL
    )

    out = list(
      opt_lambda = opt_lambda,
      opt_init = init,
      fit = run
    )

    out

  }, future.seed = TRUE)
  end_time = Sys.time()
  run_time = end_time - init_time

  names(runs) = lapply(models, function(model){
    paste0("G=", model["G"], ", M=", model["M"])
  }) %>% do.call(rbind, .)

  fit = purrr::map(runs, ~{.x$fit})
  lambda_opt = purrr::map(runs, ~{.x$opt_lambda})
  best_lambda = purrr::map(lambda_opt, ~{.x$best_lambda})
  init_opt = purrr::map(runs, ~{.x$opt_init})
  metrics = comp_metrics(fit, model_data_min = model_data_min) %>% format_metrics()

  out = list(
    fit = fit,
    init = init_opt,
    lambda_opt = lambda_opt,
    best_lambda = best_lambda,
    run_time = run_time,
    seed = seed,
    metrics = metrics,
    model_data = model_data_min
  )

  return(out)

}
