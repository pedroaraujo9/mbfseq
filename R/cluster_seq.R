#' Clustering sequence with a Bayesian functional mixture of categorical distributions
#'
#' Run the full cluster sequence modeling pipeline for varying numbers of clusters (M), including model data preparation, (optional) lambda calibration, robust initialization, model fitting, and metrics/summary extraction. Can be parallelized across cluster values.
#'
#' @param M Integer vector. The number(s) of clusters to fit for the secondary clustering w. Model will be run for each value of M.
#' @param z Integer vector. Primary cluster assignments for each subject or observation.
#' @param id Vector of subject or observation identifiers.
#' @param time Numeric vector indicating the time for each observation.
#' @param w (Optional) Vector of initial cluster assignments for the secondary clustering. If NULL, will be estimated.
#' @param lambda (Optional) Regularization or penalty hyperparameter. If NULL, will be tuned internally via grid search for each M.
#' @param n_basis Integer. Number of basis functions for time expansion (default: 10).
#' @param fixed_sd Numeric. Standard deviation for proposal or prior (default: 10).
#' @param iters Integer. Number of Markov chain or EM iterations for each fit (default: 1000).
#' @param thin Integer. Thinning of chain (default: 1; keep all samples).
#' @param burn_in Integer. Number of iterations to discard as burn-in (default: iters/2).
#' @param epsilon Numeric. Small value used as prior/jitter (default: 1).
#' @param bounds Numeric vector of length 2. Min and max lambda explored in grid search (default: c(0.01, 5)).
#' @param max_iter_start Integer. Number of steps per short chain for robust initialization (default: 30).
#' @param n_start Integer. Number of initialization chains (default: 30).
#' @param n_cores Integer. Number of CPU cores to use for parallelization across M (default: 1).
#' @param n_cores_init Integer. Number of CPU cores used for parallelization during initialization step (default: 1).
#' @param intercept Logical. Whether to include an intercept in the model (default: FALSE).
#' @param verbose Logical. If TRUE, prints progress messages (default: FALSE).
#' @param seed (Optional) Random seed for reproducibility.
#'
#' @return A list with:
#'   \describe{
#'     \item{fits}{List of fitted model outputs from clustseq_single_run, one for each value of M.}
#'     \item{metrics}{Sequence clustering evaluation metrics computed from all fits (see comp_seq_metrics).}
#'     \item{opt_lambda}{Numeric vector of optimal lambda values selected for each M.}
#'     \item{run_time}{Total run time for the procedure (as a difftime object).}
#'   }
#'
#' @details
#' For each value of M:
#'   \itemize{
#'     \item Creates model data structure via create_model_data
#'     \item If lambda is NULL, chooses an optimal lambda via calibrate_lambda (grid search)
#'     \item Finds robust initialization via clustseq_find_init (using multiple short chains)
#'     \item Fits full model via clustseq_single_run
#'   }
#' Multi-core execution (across values of M and/or for initialization) is implemented via the future and future.apply packages. Progress is optionally displayed via the progressr package.
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map_dbl
#' @seealso create_model_data, calibrate_lambda, clustseq_find_init, clustseq_single_run, comp_seq_metrics
#' @examples
#' \dontrun{
#'   time <- rep(1:3, 10)
#'   id <- rep(1:10, each = 3)
#'   z <- rep(1:2, each = 15)
#'   set.seed(123)
#'   res <- cluster_seq(M = 2:3, z = z, id = id, time = time, n_cores = 2, n_cores_init = 2,
#'                      verbose = TRUE, n_start = 3, max_iter_start = 3, iters = 10)
#'   str(res)
#' }
#' @export
cluster_seq = function(M,
                       z,
                       id,
                       time,
                       w = NULL,
                       lambda = NULL,
                       n_basis = 10,
                       fixed_sd = 10,
                       iters = 1000,
                       thin = 1,
                       burn_in = iters/2,
                       epsilon = 1,
                       bounds = c(0.01, 5),
                       max_iter_start = 30,
                       n_start = 30,
                       n_cores = 1,
                       n_cores_init = 1,
                       intercept = FALSE,
                       verbose = FALSE,
                       seed = NULL) {


  if(!is.null(seed)) set.seed(seed)

  if(n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  }else{
    future::plan(future::sequential)
  }

  init_time = Sys.time()

  fits = progressr::with_progress({

    pr = progressr::progressor(along = M)

    future.apply::future_lapply(M, function(m){

      pr(sprintf("Running M = %d", m))

      model_data = create_model_data(
        id = id,
        time = time,
        z = z,
        w = w,
        n_basis = n_basis,
        G = max(z),
        M = m,
        intercept = intercept
      )

      if(is.null(lambda)) {

        if(verbose == TRUE) {
          cat(paste0("M = ", m, " - Finding optimal lambda \n"))
        }

        opt_lambda = calibrate_lambda(
          model_data = model_data,
          bounds = bounds,
          fixed_sd = fixed_sd,
          method = "grid"
        )

        best_lambda = opt_lambda$best_lambda

      }else{

        opt_lambda = NULL
        best_lambda = lambda

      }


      if(verbose == TRUE) {
        cat(paste0("M = ", m, " - Finding init values \n"))
      }

      find_init = clustseq_find_init(
        model_data = model_data,
        n_start = n_start,
        max_iter_start = max_iter_start,
        n_cores_init = n_cores_init,
        w = w,
        lambda = best_lambda,
        n_basis = n_basis,
        fixed_sd = fixed_sd,
        epsilon = epsilon,
        intercept = intercept,
        init_list = NULL,
        verbose = FALSE,
        seed = NULL
      )

      if(verbose == TRUE) {
        cat(paste0("M = ", m, " - Fitting model \n"))
      }

      out = clustseq_single_run(
        model_data = model_data,
        w = w,
        iters = iters,
        thin = thin,
        burn_in = burn_in,
        lambda = best_lambda,
        n_basis = n_basis,
        fixed_sd = fixed_sd,
        epsilon = epsilon,
        intercept = intercept,
        init_list = find_init$init_list,
        verbose = verbose,
        seed = NULL
      )

      out$opt_lambda = opt_lambda
      out

    }, future.seed = TRUE)
  })

  future::plan(future::sequential)

  end_time = Sys.time()
  run_time = end_time - init_time

  out = list(fits = fits)
  names(out$fits) = paste0("M", M)
  out$metrics = comp_seq_metrics(out$fits)
  out$opt_lambda = purrr::map_dbl(out$fits, ~{.x$sample_list$lambda[1,1,1]})
  out$run_time = run_time

  return(out)
}
