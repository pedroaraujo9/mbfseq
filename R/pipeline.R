#' Model Fitting Pipeline for a Single G and M Combination
#'
#' Internal: Orchestrates the complete model fitting process for a single combination
#' of G (primary clusters) and M (secondary clusters). This includes creating
#' model data, calibrating lambda (if needed), finding good initial values, and
#' running the main MCMC chain.
#'
#' @param cluster_dim Numeric vector of length 2, named "G" and "M", specifying
#'   the number of primary and secondary clusters for this run.
#' @param z Optional integer vector. Primary cluster assignments for each subject or observation.
#' @param w Optional integer vector. Secondary cluster assignments for each subject or observation.
#' @param x Optional matrix or data.frame of covariates (columns = covariates, rows = observations).
#' @param id Vector of subject or group identifiers.
#' @param time Numeric vector indicating the time for each observation.
#' @param iters Integer. Total number of MCMC iterations for the final run.
#' @param burn_in Integer. Number of burn-in iterations for the final run.
#' @param thin Integer. Thinning interval for MCMC sampling in the final run.
#' @param lambda Optional numeric value for the regularization parameter. If `NULL`,
#'   the optimal lambda is estimated via `calibrate_lambda`.
#' @param n_basis Integer. Number of spline basis functions (includes intercept).
#' @param init_list Optional list of initial values for model parameters. If `NULL`,
#'   initial values are found using `find_init`.
#' @param config List of additional configuration parameters, including options
#'   for lambda calibration and initialization search (see `calibrate_lambda` and `find_init`).
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param seed Optional integer seed for reproducibility. If `NULL`, a random seed is generated.
#'
#' @return A list containing:
#'   \describe{
#'     \item{opt_lambda}{Result from lambda calibration (see `calibrate_lambda`), or a list
#'                       with the provided lambda.}
#'     \item{opt_init}{Result from the initialization search (see `find_init`).}
#'     \item{fit}{Result from the final model run (see `single_run`).}
#'     \item{seed}{The random seed used for this pipeline run.}
#'   }
#'
#' @details
#' This function serves as an internal helper for `mbfseq_fit`, handling the sequential
#' steps for a single (G, M) model configuration. It first calls `create_model_data`
#' to prepare the data structures. If `lambda` is not provided, it performs lambda
#' tuning using `calibrate_lambda`. Then, it finds robust initial values using
#' `find_init`. Finally, it executes the main MCMC chain with `single_run` using
#' the found or provided initial values and lambda. Progress messages are printed
#' if `verbose` is `TRUE`.
#'
#' @seealso \code{\link{mbfseq_fit}}, \code{\link{create_model_data}}, \code{\link{calibrate_lambda}}, \code{\link{find_init}}, \code{\link{single_run}}
#' @keywords internal
pipeline = function(cluster_dim,
                    z,
                    w,
                    x,
                    id,
                    time,
                    iters,
                    burn_in,
                    thin,
                    lambda,
                    n_basis,
                    init_list,
                    config = list(
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
                    seed = NULL) {

  if(is.null(seed)) seed = sample(1:10000, size = 1)
  set.seed(seed)

  g = cluster_dim["G"]
  m = cluster_dim["M"]

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

  # find optimal lambda if lambda is not provided
  if(is.null(lambda)) {

    if(verbose) cat(paste0("G = ", g, ", M = ", m, " - Finding lambda\n"))

    opt_lambda = calibrate_lambda(
      z = z,
      w = w,
      model_data = model_data,
      config = config
    )

  }else{

    opt_lambda = list(best_lambda = lambda)

  }

  # find inits
  if(verbose) cat(paste0("G = ", g, ", M = ", m, " - Finding inits\n"))

  init = find_init(
    n_start = config$n_start,
    iters = config$n_start_iters,
    n_cores = config$n_start_cores,
    model_data = model_data,
    lambda = opt_lambda$best_lambda,
    init_list = init_list,
    priors = config,
    seed = NULL
  )

  # final run
  run = single_run(
    model_data = model_data,
    lambda = opt_lambda$best_lambda,
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
    fit = run,
    seed = seed
  )

  return(out)

}
