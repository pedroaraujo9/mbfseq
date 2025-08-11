#' Find Initial Values via Multiple Restarts
#'
#' Performs multiple runs of the model with different initializations to find
#' good starting values for the main MCMC chain. Runs are executed in parallel
#' when possible.
#'
#' @param n_start Integer specifying number of initialization runs (default: 30)
#' @param iters Integer specifying number of iterations per run (default: 50)
#' @param n_cores Integer specifying number of cores for parallel execution (default: 1)
#' @param model_data List containing all required model data components
#' @param lambda Numeric penalty parameter for regularization (default: 1)
#' @param init_list Optional list of initial values for parameters. If NULL,
#'                  random initialization is used (default: NULL)
#' @param priors List containing prior parameters:
#' \describe{
#'   \item{epsilon_w}{Precision for w prior (default: 1)}
#'   \item{beta_sd}{SD for non-penalized coefficients (default: sqrt(10))}
#'   \item{mu_sd}{SD for cluster means prior (default: sqrt(10))}
#'   \item{sigma_a}{Shape parameter for sigma prior (default: 1)}
#'   \item{sigma_b}{Rate parameter for sigma prior (default: 1)}
#' }
#' @param seed Optional random seed for reproducibility (default: NULL)
#'
#' @return A list containing:
#' \describe{
#'   \item{logpost}{Vector of mean log-posterior values for each run}
#'   \item{init_list}{List of parameter values from the best run}
#'   \item{run_time}{Total execution time for all runs}
#' }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Runs `n_start` independent chains in parallel (if `n_cores > 1`)
#'   \item For each chain:
#'   \itemize{
#'     \item Runs for `iters` iterations with burn-in of `iters/2`
#'     \item Thins samples by factor of 2
#'   }
#'   \item Selects the chain with highest mean log-posterior value
#'   \item Returns the final parameter values from the best chain as initial values
#' }
#'
#' Parallel execution uses `future` and `future.apply` packages.
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom purrr map_dbl map
#' @importFrom stats runif
#'
#' @examples
#' \dontrun{
#' # Example usage (not run)
#' # Find initial values with 10 restarts using 2 cores
#' init_results <- find_init(
#'   n_start = 10,
#'   iters = 100,
#'   n_cores = 2,
#'   model_data = data,
#'   lambda = 1.5,
#'   priors = priors,
#'   seed = 42
#' )
#' }
#'
#' @seealso \code{\link{single_run}} for the main MCMC function that uses these initial values
#' @export
find_init = function(n_start = 30,
                     iters = 30,
                     n_cores = 1,
                     model_data,
                     lambda = 1,
                     init_list = NULL,
                     priors = list(
                       epsilon_w = 1,
                       beta_sd = sqrt(10),
                       mu_sd = sqrt(10),
                       sigma_a = 1,
                       sigma_b = 1
                     ),
                     seed = NULL) {

  if(is.null(seed)) seed = sample(1:10000, size = 1)
  set.seed(seed)

  ## define session
  if(n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  }else{
    future::plan(future::sequential)
  }

  init_time = Sys.time()
  runs = future.apply::future_lapply(1:n_start, function(i){

    single_run(
      model_data = model_data,
      lambda = lambda,
      iters = iters,
      burn_in = iters/2,
      thin = 2,
      init_list = init_list,
      priors = priors,
      verbose = FALSE,
      seed = NULL
    )

  }, future.seed = TRUE)
  end_time = Sys.time()

  logpost = map_dbl(runs, ~{mean(.x$logpost[, "logpost"])})
  best_run = runs[[which.max(logpost)]]
  init_list = purrr::map(best_run$sample_list, ~{.x |> tail(1)})

  out = list(
    logpost = logpost,
    init_list = init_list,
    run_time = end_time - init_time
  )

  return(out)

}
