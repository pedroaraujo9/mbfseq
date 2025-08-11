#' Update Lambda Parameters
#'
#' Internal: Update the lambda (penalty) parameters for each group and class using a Gamma or inverse-Gamma prior, based on the current values of the alpha coefficients.
#'
#' @param alpha Numeric matrix of current alpha coefficients (basis coefficients), with rows corresponding to basis functions and columns to groups/classes.
#' @param a Numeric. Shape parameter for the Gamma prior (default: 1).
#' @param b Numeric. Rate parameter for the Gamma prior (default: 1).
#' @param model_data List. Model data structure (from create_model_data), must include G, M, n_basis, notpen_index, and other relevant fields.
#' @param grouped_lambda Logical. If TRUE, use grouped lambda update (default: FALSE).
#' @param CUST Logical. Custom mode (default: FALSE).
#'
#' @return Matrix of updated lambda values with dimensions M (groups) by G-1 (classes).
#'
#' @details
#' For each class (except the last), this function aggregates the squared alpha coefficients (excluding non-penalized indices) by group, and samples new lambda values from the (inverse-)Gamma posterior. The update is performed for each group and class. If grouped_lambda is TRUE, a grouped update is performed. For penalized splines (basis_type == "ps"), a different update is used based on the S matrix.
#'
#' @importFrom magrittr %>%
#' @importFrom stats aggregate rgamma
#' @importFrom extraDistr rinvgamma
#' @keywords internal
update_lambda = function(alpha,
                         a = 1,
                         b = 1,
                         model_data,
                         grouped_lambda = FALSE,
                         CUST = FALSE) {
  G = model_data$G
  M = model_data$M
  n_basis = model_data$n_basis
  notpen_index = model_data$notpen_index
  order = model_data$order
  id_effect = model_data$include_id_effect
  n_id = model_data$n_id
  basis_type = model_data$basis_type
  newD = model_data$newD

  if(basis_type != "ps") {

    if(id_effect == TRUE) {
      alpha = alpha[-c(1:n_id), ]
    }


    if(grouped_lambda == FALSE) {

      #notpen_index = gen_notpen_index(n_basis = n_basis, M = M, order = 2)
      #n_params = rep(n_basis - 2, M)
      n_params = rep(n_basis - order, M)

      new_lambda = lapply(1:(G-1), function(g){

        #sq_sum = (alpha[-notpen_index, g]^2) |> matrix(ncol = M) |> colSums()

        alpham = alpha[-notpen_index, g] |> matrix(ncol = M)
        sq_sum = diag(t(alpham) %*% newD[-1, ][, -1] %*% alpham)

        extraDistr::rinvgamma(n = M, alpha = a + n_params/2, beta = b + sq_sum/2)

      }) %>% do.call(cbind, .)

    }else{

      n_params = rep((G-1)*(n_basis - order), M)

      sq_sum = lapply(1:(G-1), function(g){

        aggregate(
          alpha[-notpen_index, g] ~ rep(1:M, each = n_basis-order),
          FUN = function(x) sum(x^2)
        )[, -1]

      }) %>%
        do.call(cbind, .) %>%
        rowSums()

      new_lambda = 1/stats::rgamma(n = M, a + n_params, b + sq_sum/2) |>
        rep(G-1) |>
        matrix(nrow = M)
    }
  }else{

    S = model_data$S
    n_params = rep(n_basis, M)
    new_lambda = lapply(1:(G-1), function(g){

      alpha_m = matrix(alpha[, g], nrow = n_basis)
      sq_sum = diag(t(alpha_m) %*% S %*% alpha_m)

      1/stats::rgamma(n = M, a + n_params/2, b + sq_sum/2)

    }) %>% do.call(cbind, .)

  }

  return(new_lambda)

}
