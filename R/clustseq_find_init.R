#' Find Robust Initial Values for Cluster Sequence Model
#'
#' Run multiple short chains for the cluster sequence model in parallel (or sequentially) to search for robust initial parameter values and cluster assignments. Each short chain fits using clustseq_single_run; the run with the highest mean log posterior is selected, and its last sample is used for initialization.
#'
#' @param model_data A model_data object from create_model_data, with all basis, cluster, and covariate fields needed for model fitting.
#' @param n_start Integer. Number of short MCMC/EM chains to run for initialization (default: 30).
#' @param max_iter_start Integer. Number of iterations in each short initialization chain (default: 30).
#' @param n_cores_init Integer. Number of CPU cores for parallelization across initial chains (default: 1).
#' @param w (Optional) Initial secondary cluster assignment vector for all initialization chains.
#' @param lambda (Optional) Penalty or regularization hyperparameter to use during initialization.
#' @param n_basis Integer. Number of basis functions for modeling time (default: 10).
#' @param fixed_sd Numeric. Fixed standard deviation for proposals/priors (default: 10).
#' @param epsilon Numeric. Small value for priors/jitter (default: 1).
#' @param intercept Logical. Whether to fit an intercept term in the model (default: FALSE).
#' @param init_list (Optional) List of specific initial parameter values for the chains.
#' @param verbose Logical. If TRUE, prints progress for each short chain (default: FALSE).
#' @param seed (Optional) Integer random seed for reproducible initialization.
#'
#' @return A list with:
#'   \describe{
#'     \item{best_fit}{The full output of clustseq_single_run from the initial chain with the highest mean log posterior.}
#'     \item{init_list}{List of final parameter values (the last sample from best_fit$sample_list) for warm starting further fitting.}
#'     \item{logpost_mean}{Numeric vector with the mean log posterior across all initialization chains.}
#'   }
#'
#' @details
#' Uses the future and future.apply packages for multicore parallelization across initial chains. Especially suited for mixture/cluster models where initialization is sensitive and local optima are a risk. The function always restores the R session's future plan to sequential on exit. This function is intended for internal or parallel use and is not exported to end users.
#'
#' @importFrom future plan multisession sequential reset
#' @importFrom future.apply future_lapply
#' @importFrom purrr map_dbl map
#' @examples
#' \dontrun{
#'   # Example: Robust initialization for cluster sequence modeling
#'   set.seed(2023)
#'   time <- rep(1:3, 10)
#'   id <- rep(1:10, each = 3)
#'   z <- rep(1:2, each=15)
#'   w <- rep(1:2, length.out=30)
#'   model_data <- create_model_data(time = time, id = id, z = z, w = w)
#'
#'   # Parallel, short runs for initialization
#'   init <- clustseq_find_init(
#'     model_data = model_data,
#'     n_start = 5,
#'     max_iter_start = 5,
#'     n_cores_init = 2
#'   )
#'   str(init)
#' }
clustseq_find_init = function(model_data,
                              n_start = 30,
                              max_iter_start = 30,
                              n_cores_init = 1,
                              w = NULL,
                              lambda = NULL,
                              n_basis = 10,
                              fixed_sd = 10,
                              epsilon = 1,
                              intercept = FALSE,
                              init_list = NULL,
                              verbose = FALSE,
                              seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  if(n_cores_init > 1) {
    future::plan(future::multisession, workers = n_cores_init)
  }else{
    future::plan(future::sequential)
  }

  fits = future.apply::future_lapply(1:n_start, function(i){

    clustseq_single_run(
      model_data = model_data,
      w = w,
      iters = max_iter_start,
      thin = 1,
      burn_in = 1,
      lambda = lambda,
      n_basis = n_basis,
      fixed_sd = fixed_sd,
      epsilon = epsilon,
      intercept = intercept,
      init_list = init_list,
      verbose = verbose
    )

  }, future.seed = TRUE)

  if(n_cores_init > 1) {

    future::plan(future::multisession, workers = n_cores_init)

    fits = future.apply::future_lapply(1:n_start, function(i){

      clustseq_single_run(
        model_data = model_data,
        w = w,
        iters = max_iter_start,
        thin = 1,
        burn_in = 1,
        lambda = lambda,
        n_basis = n_basis,
        fixed_sd = fixed_sd,
        epsilon = epsilon,
        intercept = intercept,
        init_list = init_list,
        verbose = verbose
      )

    }, future.seed = TRUE)
  }

  future::reset(future::multisession)
  future::plan(future::sequential)

  logpost_mean = purrr::map_dbl(fits, ~{
    mean(.x$sample_list$w_logpost)
  })

  best_fit = fits[[which.max(logpost_mean)]]
  init_list = purrr::map(best_fit$sample_list, ~{.x |> tail(1)})

  out = list(
    best_fit = best_fit,
    init_list = init_list,
    logpost_mean = logpost_mean
  )

  return(out)

}
