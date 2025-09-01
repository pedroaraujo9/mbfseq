#' Find MAP Estimates via Hierarchical Clustering and Stan Optimization
#'
#' Computes a maximum a posteriori estimate (MAP) of model parameters and cluster structure using hierarchical clustering for initialization and Stan optimization for parameter estimation.
#'
#' @param z Sequence vector for clustering initialization.
#' @param w Cluster assignments for initialization.
#' @param lambda Penalty or regularization hyperparameter(s) passed to Stan and the optimizer.
#' @param stan_model A compiled Stan model object (see `rstan`).
#' @param model_data A `model_data` object (from `create_model_data`) containing all basis, grouping, and index info.
#' @param fixed_sd Numeric. Fixed standard deviation for prior/penalty in the model. Default is 10.
#'
#' @return A list with components:
#'   - `opt_fit`: The result of Stan MAP optimization (see `rstan::optimizing`).
#'
#' @details
#' Uses cluster assignments and basis expansions to generate the Stan data list. Runs MAP optimization using the LBFGS algorithm via `rstan`. Intended for internal use.
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

  M = length(unique(w))
  n_basis = model_data$n_basis

  w = w |> factor(levels = 1:M)
  if(M > 1) {
    W = model.matrix(~ - 1 + w)
  }else{
    W = matrix(1, nrow = length(w), ncol = 1)
  }


  D = model_data$nD %>% diag() %>% .[-1]

  if(M == 1) {
    notpen_index = array(1)
  }else{
    notpen_index = model_data$notpen_index
  }

  data_list = list(
    n = model_data$n,
    G = model_data$G,
    M = M,
    n_basis = model_data$n_basis,
    order = model_data$order,
    notpen_index = notpen_index,
    spline_index = cbind(matrix(1:((M * n_basis)), ncol = M)[-1, ]),
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
    opt_fit = opt_fit,
    data_list = data_list
  )

  return(out)

}

#' Calibrate Lambda for Penalization via Grid Search
#'
#' Optimizes the penalty hyperparameter lambda for the cluster sequence model by grid search, using held-out likelihood as the objective.
#'
#' @param z Optional. Sequence vector for clustering initialization. If `NULL`, will use `find_init` to generate initial `z` and `w`.
#' @param w Optional. Cluster assignments for initialization. If `NULL`, will be computed from `z` or `find_init`.
#' @param model_data A `model_data` object as returned by `create_model_data`, with basis, clustering, and index fields.
#' @param single_group Logical. If `TRUE`, assumes a single group for clustering. Default is `FALSE`.
#' @param config List of configuration options:
#'   - `bounds`: Numeric vector of length 2. Range of lambda to explore. Default is `c(0.01, 10)`.
#'   - `n_points`: Number of grid points. Default is 20.
#'   - `n_start_iters`: Number of initialization iterations for `find_init`. Default is 20.
#'   - `lambda_start`: Initial lambda value for `find_init`. Default is 1.
#'   - `epsilon_w`, `beta_sd`, `mu_sd`, `sigma_a`, `sigma_b`: Priors and scaling constants.
#'
#' @return A list with the following elements:
#'   - `fit`: A data frame with lambda values and penalized log-likelihood scores.
#'   - `best_lambda`: The optimal value of lambda found.
#'   - `w`: Cluster assignments used for optimization.
#'   - `z`: Sequence vector used for optimization.
#'   - `lambda_plot`: `ggplot2` object for visualizing penalized likelihood versus lambda.
#'
#' @details
#' Loads and reuses a Stan model template from the package's `extdata` folder. Uses `find_map` to evaluate the penalized likelihood for each value of lambda. Intended for internal use.
#'
#' @importFrom purrr map_dbl
#' @importFrom ggplot2 ggplot aes geom_point geom_line
#' @keywords internal
calibrate_lambda = function(z = NULL,
                            w = NULL,
                            model_data,
                            single_group = FALSE,
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

    if(single_group == TRUE) {
      w = rep(1, model_data$n_id)
    }else{
      w = hclust(z_ham_dist, method = "ward.D") |> cutree(k = model_data$M)
    }

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

    if(is.null(z)) z = run$init_list$z[1,]
    if(is.null(w)) w = run$init_list$w[1,]
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

