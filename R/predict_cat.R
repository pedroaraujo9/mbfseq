#' Predict categorical outcomes from model coefficients
#'
#' Internal function to generate predictions from a fitted categorical regression model,
#' returning either the linear predictor, full probability vector, or probability of the observed category.
#'
#' @param alpha A numeric matrix of model coefficients.
#' @param z A vector of observed categorical outcomes (integer-encoded, 1 to G).
#' @param w_vec A named vector with the cluster allocation for each observation.
#' @param time_seq A vector of time indices corresponding to each observation.
#' @param model_data A list containing model structure elements.
#' @param type Character string indicating the type of prediction. Options are \code{"linear"}, \code{"probability"}, or \code{"category probability"}.
#' @param intercept Logical; whether to include an intercept in the design matrix.
#'
#' @return A numeric matrix or vector:
#'   \itemize{
#'     \item If \code{type = "linear"}, returns the linear predictor.
#'     \item If \code{type = "probability"}, returns the softmax-transformed probabilities for each category.
#'     \item If \code{type = "category probability"}, returns the probability associated with the observed category for each observation.
#'   }
#'
#' @importFrom mclust softmax
#' @keywords internal
#' @noRd
predict_cat = function(alpha,
                       z,
                       w_vec,
                       time_seq,
                       model_data,
                       type = "probability",
                       intercept = FALSE) {

  G = model_data$G
  M = model_data$M
  n_basis = model_data$n_basis
  B_unique = model_data$B_unique
  id_unique = model_data$id_unique

  B = B_unique[time_seq, , drop = FALSE]
  B_expand = B[, rep(1:n_basis, times = M)]

  X = gen_reg_matrix(
    w_vec = w_vec,
    M = M,
    B_expand = B_expand,
    intercept = intercept
  )

  linear_pred = X %*% alpha
  theta = mclust::softmax(linear_pred)

  if(type == "linear") {

    return(linear_pred)

  }else if(type == "probability") {

    return(theta)

  }else if(type == "category probability") {

    Z = gen_dummy(z, n_cat = G, intercept = FALSE)
    cat_theta = rowSums(Z * theta)
    return(cat_theta)

  }

}
