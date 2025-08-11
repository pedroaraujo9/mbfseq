gv = c(
  c(".", "nvars", "w", "time", "li", "ui", "id", "logp_z", "lambda")
)

utils::globalVariables(gv)

#' Generate Dummy (Indicator) Matrix
#'
#' Internal: Create a dummy variable (indicator) matrix for a categorical variable, with optional intercept.
#'
#' @param x Integer vector of categories (1-based).
#' @param n_cat Integer. Number of categories.
#' @param intercept Logical. Include intercept column? Default is FALSE.
#'
#' @return Numeric matrix of dummy variables.
#' @export
gen_dummy = function(x, n_cat, intercept = FALSE) {

  x = factor(x, levels = 1:n_cat)

  if(intercept == FALSE) {

    if(n_cat == 1) {
      X = cbind(rep(1, length(x)))

    }else{
      X = model.matrix( ~ -1 + x)
      X = X[1:nrow(X), ]
    }

  }else{

    if(n_cat == 1) {
      X = cbind(rep(1, length(x)))

    }else{
      X = model.matrix( ~ x)
      X = X[1:nrow(X), ]
    }

  }

  colnames(X) = NULL
  rownames(X) = NULL

  return(X)
}

#' Indices of Non-Penalized Basis Functions
#'
#' Internal: Generate indices for non-penalized basis functions across groups/classes.
#'
#' @param n_basis Integer. Number of basis functions per group.
#' @param M Integer. Number of groups/classes.
#' @param order Integer. Order of difference for penalization.
#'
#' @return Integer vector of indices.
#' @importFrom magrittr %>%
#' @keywords internal
gen_notpen_index = function(n_basis, M, order = 1) {
  notpen_index = lapply(1:M, function(k){
    (k -1) * n_basis + 1:(order)
  }) %>% do.call(c, .)
  return(notpen_index)
}

#' Indices for Grouped Basis Functions
#'
#' Internal: Create a list of index vectors for basis functions by group.
#'
#' @param n_basis Integer. Number of basis functions per group.
#' @param M Integer. Number of groups.
#'
#' @return List of integer vectors, one per group.
#' @keywords internal
gen_basis_index = function(n_basis, M) {
  lapply(1:M, function(k){
    (((k-1)*(n_basis) + 1):(k * (n_basis)))
  })
}

#' Generate Inverse Covariance Matrices for Penalized Coefficients
#'
#' Internal: Create a list of diagonal inverse covariance matrices for each class (except last).
#'
#' @param lambda Numeric matrix or vector of penalty parameters.
#' @param model_data List. Model structure (see create_model_data).
#' @param fixed_sd Numeric. SD for non-penalized coefficients (default 100).
#'
#' @return List of diagonal matrices.
#' @keywords internal
gen_inv_cov = function(lambda, model_data, fixed_sd = 100) {

  G = model_data$G
  M = model_data$M
  n_basis = model_data$n_basis
  notpen_index = model_data$notpen_index
  n_id = model_data$n_id
  nD = model_data$nD
  lambda = cbind(lambda)

  var_vec = rep((lambda), each = n_basis * M)
  D = rep(1/diag(nD), times = M)
  var_vec = var_vec * D
  var_vec[notpen_index] = (fixed_sd^2)
  prec_matrix = diag(1/var_vec)

  inv_cov_list = lapply(1:(G-1), function(g){prec_matrix})
  return(inv_cov_list)

}

#' Regularized Design Matrix for Cluster-Specific Effects
#'
#' Internal: Construct design matrix for cluster-specific basis coefficients.
#'
#' @param w_vec Integer vector. Cluster assignment for each row.
#' @param M Integer. Number of clusters.
#' @param B_expand Matrix. Basis expansion for all clusters.
#' @param intercept Logical. Include intercept? Default is FALSE.
#'
#' @return Numeric matrix, same shape as B_expand.
#' @keywords internal
gen_reg_matrix = function(w_vec, M, B_expand, intercept = FALSE) {

  n_basis = ncol(B_expand)/(M)

  W = gen_dummy(w_vec, n_cat = M, intercept = intercept)
  W_expand = W[, rep(1:M, each = n_basis)]

  X = B_expand * W_expand

  return(X)
}

#' Posterior Summary Statistic by Group
#'
#' Internal: Compute summary statistic (e.g., mean) across posterior samples for each group.
#'
#' @param sample 3D numeric array: iterations x items x groups.
#' @param stat_function Function to apply (e.g., mean).
#' @param ... Additional arguments to stat_function.
#'
#' @return Numeric matrix: items x groups.
#' @importFrom magrittr %>%
#' @keywords internal
comp_post_stat = function(sample, stat_function = mean, ...) {
  est = lapply(1:dim(sample)[3], function(g){
    sample[, , g] |> apply(MARGIN = 2, FUN = stat_function, ...)
  }) %>%
    do.call(cbind, .)
  return(est)
}

#' Most Frequent Class Assignment (Posterior Mode)
#'
#' Internal: Compute the most frequent class for each item across samples.
#'
#' @param sample Matrix or array: iterations x items.
#'
#' @return Integer vector of modal class assignments.
#' @keywords internal
comp_class = function(sample) {

  class = apply(
    sample,
    MARGIN = 2,
    FUN = function(x) x |> table() |> which.max() |> names() |> as.integer()
  )

  return(class)

}

#' Convert Posterior Summary Matrix to Long Format
#'
#' Internal: Convert summary matrix to tidy long-format tibble for plotting/analysis.
#'
#' @param summ_matrix Numeric matrix. Posterior summary (e.g., means).
#' @param col_name Character. Name for summary column.
#' @param time Numeric/integer vector. Time points.
#' @param w Integer/factor vector. Group labels.
#'
#' @return Tibble in long format: time, w, cat, summary.
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom stringr str_replace
#' @keywords internal
trans_summ_prob = function(summ_matrix, col_name, time, w) {
  summ_matrix |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::mutate(time = time, w = w) |>
    tidyr::gather(cat, {{col_name}}, -time, -w) |>
    dplyr::mutate(cat = str_replace(cat, "V", ""))
}

#' Compute AICM and BICM for Sequential Clustering Fits
#'
#' Internal: Compute AICM, BICM, and log-posterior stats for a list of fits.
#'
#' @param fits List of fitted model objects.
#'
#' @return List with numeric vectors: AICM, BICM, logp_var, logp_avg.
#' @importFrom stats var
#' @importFrom purrr map_dbl
#' @keywords internal
comp_seq_metrics = function(fits) {

  metrics = lapply(fits, FUN = function(fit){

    n = fit$model_data$n_id
    logp = fit$sample_list$w_logpost

    metrics = list(
      AICM = -2*(mean(logp) - stats::var(logp)),
      BICM = -2*(mean(logp) - (log(n) - 1)*stats::var(logp)),
      logp_var = var(logp),
      logp_avg = mean(logp)
    )
    return(metrics)
  })

  out = list(
    AICM = purrr::map_dbl(metrics, ~{.x$AICM}),
    BICM = purrr::map_dbl(metrics, ~{.x$BICM}),
    logp_var = purrr::map_dbl(metrics, ~{.x$logp_var}),
    logp_avg = purrr::map_dbl(metrics, ~{.x$logp_avg})
  )

  return(out)
}

filter_array = function(array, burn_in, thin) {

  if(!is.null(array)) {

    array_dim = dim(array)
    dim_len = length(array_dim)
    iters = array_dim[1]
    post_iters = seq(burn_in, iters, thin)

    if(dim_len == 1) {
      post_array = array[post_iters]
    }else if(dim_len == 2) {
      post_array = array[post_iters, ]
    }else if(dim_len == 3) {
      post_array = array[post_iters,,]
    }

  }else{
    post_array = NULL
  }

  return(post_array)
}

gen_sample_array = function(iters, dimension, sampler = NULL, init = NULL) {

  dim_len = length(dimension)
  sample_array = array(dim = c(iters, dimension))


  if(is.null(init)) {
    init = sampler(prod(dimension))
  }

  if(dim_len == 0) {
    sample_array[1] = init
  }else if(dim_len == 1) {
    sample_array[1,] = init
  }else if(dim_len == 2) {
    sample_array[1,,] = init
  }

  return(sample_array)

}

