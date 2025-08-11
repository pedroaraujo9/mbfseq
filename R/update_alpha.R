#' Update Alpha Coefficient Matrix
#'
#' Internal: Update the matrix of coefficients (alpha) for the clustering model using current assignments and regularization/penalization parameters. Used in the MCMC or EM loop to sample or optimize new alpha values for each group.
#'
#' @param alpha Matrix. Current values of the alpha coefficients (basis × groups).
#' @param lambda Regularization or penalty hyperparameter(s) for the prior or penalty (format as required by gen_inv_cov).
#' @param z Integer vector. Primary cluster assignment for each subject or observation.
#' @param w Integer vector. Secondary cluster assignment for each subject or observation.
#' @param model_data List. Model data structure (from create_model_data), must include all required matrices, indices, and group counts.
#' @param fixed_sd Numeric. Fixed standard deviation for the prior or penalty in the update step (default: 10).
#' @param intercept Logical. Whether to include an intercept in the model (default: FALSE).
#'
#' @return Updated matrix of alpha coefficients (basis × groups).
#'
#' @details
#' This function cycles through G-1 groups (leaving one as baseline/reference), updating coefficients for each group. It constructs dummy matrices for group assignments (Z and W), forms the block design matrix via kronecker product, and calls internal helpers (gen_inv_cov, gen_dummy, augment_alpha) to perform the update. The update is typically used within a Gibbs sampler or EM algorithm for mixture or clustering models with basis expansions.
#'
#' @seealso clustseq_single_run, create_model_data
#' @keywords internal
update_alpha = function(alpha,
                        lambda,
                        z,
                        w,
                        model_data,
                        fixed_sd = 10,
                        intercept = FALSE) {

  G = model_data$G
  M = model_data$M
  B_expand = model_data$B_expand
  id = model_data$id_time_df$id
  basis_index = model_data$basis_index
  n_basis = model_data$n_basis
  B_unique = model_data$B_unique
  notpen_index = model_data$notpen_index

  alpha = rbind(alpha)

  inv_cov_list = gen_inv_cov(
    lambda,
    model_data = model_data,
    fixed_sd = fixed_sd
  )

  W = gen_dummy(w, M, intercept = intercept)
  X = kronecker(W, B_unique)
  Z = gen_dummy(z, n_cat = G, intercept = FALSE)

  for(g in 1:(G-1)) {

    linear_pred = X %*% alpha

    alpha[, g] = augment_alpha(
      g = g,
      Z = Z,
      X = X,
      linear_pred = linear_pred,
      inv_cov_list = inv_cov_list
    )
  }

  return(alpha)

}
