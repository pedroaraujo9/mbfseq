#' Find MAP Estimates via Hierarchical Clustering and Stan Optimization
#'
#' Internal: Compute a maximum a posteriori estimate (MAP) of model parameters and cluster structure using hierarchical clustering for initialization and Stan optimization for parameter estimation.
#'
#' @param lambda Penalty or regularization hyperparameter(s) passed to Stan and the optimizer.
#' @param model A compiled Stan model object (see rstan).
#' @param model_data A model_data object (from create_model_data) containing all basis, grouping, and index info.
#' @param fixed_sd Numeric. Fixed standard deviation for prior/penalty in the model (default: 10).
#'
#' @return A list with components:
#'   \describe{
#'     \item{opt_fit}{The result of Stan MAP optimization (see rstan::optimizing).}
#'     \item{w_ward}{Vector of clusters assigned by Ward hierarchical clustering on Hamming distances (z_ham_dist in model_data).}
#'   }
#'
#' @details
#' Runs hierarchical clustering (hclust with Ward's method) on the Hamming distance between sequences to obtain initial clusters. Clusters are used to generate a model matrix (W), which, along with basis expansions and model indices, form the Stan data list. Runs MAP optimization using the LBFGS algorithm via rstan. Intended for internal use.
#'
#' @importFrom stats hclust cutree model.matrix
#' @importFrom rstan optimizing
#' @importFrom magrittr %>%
#' @keywords internal
find_map = function(z,
                    w,
                    lambda,
                    stan_model,
                    model_data,
                    fixed_sd = 10) {

  M = model_data$M
  n_basis = model_data$n_basis

  w = w |> factor(levels = 1:M)
  W = model.matrix(~ - 1 + w)

  D = model_data$nD %>% diag() %>% .[-1]

  data_list = list(
    n = model_data$n,
    G = model_data$G,
    M = M,
    n_basis = model_data$n_basis,
    order = model_data$order,
    notpen_index = model_data$notpen_index,
    spline_index = matrix(1:((M * n_basis)), ncol = M)[-1, ],
    X = kronecker(W, model_data$B_unique),
    lambda = lambda,
    z = z,
    fixed_sd = fixed_sd,
    D = D
  )

  opt_fit = rstan::optimizing(
    object = stan_model,
    data = data_list,
    algorithm = "LBFGS",
    as_vector = FALSE,
    iter = 1000
  )

  out = list(
    opt_fit = opt_fit
  )

  return(out)

}

#' Calibrate Lambda for Penalization via Bayesian or Grid Search
#'
#' Internal: Optimize the penalty hyperparameter lambda for your cluster sequence model by either Bayesian optimization or grid search, using held-out likelihood as the objective.
#'
#' @param bounds Numeric vector of length 2. Range of lambda to explore (default: c(0.01, 10)).
#' @param model_data A model_data object as returned by create_model_data, with basis, clustering, and index fields.
#' @param fixed_sd Numeric. Fixed standard deviation for priors/penalties in the Stan model (default: 10).
#' @param method Character. "bayes" for Bayesian optimization via the ParBayesianOptimization package, or "grid" for a regular grid search.
#' @param options List. Options passed to the optimizer.
#'   For the "bayes" method:
#'     - init_points: Number of random initialization points (default: 5)
#'     - iters_n: Number of Bayesian optimization iterations (default: 10)
#'     - acq: Acquisition function (default: "ucb")
#'     - kappa: Exploration factor for the acquisition function (default: 2.57)
#'   For the "grid" method:
#'     - length_out: Number of grid points (default: 20)
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{fit}{The object returned by the optimizer (a bayesOpt object for the "bayes" method, or a data.frame for the "grid" method).}
#'     \item{best_lambda}{The optimal value of lambda found.}
#'     \item{lambda_plot}{(Grid method only) ggplot2 object for visualizing penalized likelihood versus lambda.}
#'   }
#'
#' @details
#' Loads and reuses a Stan model template from the package's extdata folder. Uses find_map to evaluate the penalized likelihood for each value of lambda. Bayesian optimization is powered by the ParBayesianOptimization package. For grid search, the output includes a plot of penalized score versus lambda (see lambda_plot). Intended for internal use.
#'
#' @importFrom ParBayesianOptimization bayesOpt getBestPars
#' @importFrom purrr map_dbl
#' @importFrom ggplot2 ggplot aes geom_point geom_line
#' @keywords internal
calibrate_lambda = function(z = NULL,
                            w = NULL,
                            model_data,
                            config = list(
                              bounds = c(0.01, 10),
                              n_points = 20,
                              n_start_iters = 20,
                              lambda_start = 1,
                              epsilon_w = 1,
                              beta_sd = sqrt(10),
                              mu_sd = sqrt(10),
                              sigma_a = 1,
                              sigma_b = 1
                            )) {

  model_path = system.file(
    "extdata", "model-stan.rds", package = "mbfseq"
  )

  stan_model = readRDS(model_path)

  if(!is.null(z) & is.null(w)) {
    z_ham_dist = compute_hamming(z = z, model_data = model_data)
    w = hclust(z_ham_dist, method = "ward.D") |> cutree(k = model_data$M)

  }else if(is.null(z)){

    run = find_init(
      n_start = config$n_points,
      iters = config$n_start_iter,
      n_cores = 1,
      model_data = model_data,
      lambda = config$lambda_start,
      init_list = NULL,
      priors = config,
      seed = NULL
    )

    z = run$init_list$z[1,]
    w = run$init_list$w[1,]

  }

  lambda_grid = seq(config$bounds[1], config$bounds[2], length.out = config$n_points)

  obj_fun = function(lambda) {
    fit = find_map(
      z = z,
      w = w,
      lambda = lambda,
      stan_model = stan_model,
      model_data = model_data,
      fixed_sd = config$beta_sd
    )

    return(list(Score = fit$opt_fit$par$penal_ll))

  }

  penal = purrr::map_dbl(lambda_grid, ~{
    obj_fun(.x)$Score
  })

  best_lambda = lambda_grid[which.max(penal)]

  out = list(
    fit = data.frame(lambda = lambda_grid, log_penal = penal),
    best_lambda = best_lambda,
    w = w,
    z = z
  )

  out$lambda_plot = ggplot2::ggplot(out$fit, ggplot2::aes(x=lambda, y=log_penal)) +
    ggplot2::geom_point() +
    ggplot2::geom_line()

  return(out)

}

