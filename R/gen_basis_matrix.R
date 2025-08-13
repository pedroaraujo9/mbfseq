#' Generate Penalized Spline Basis Matrix for Time Covariate
#'
#' Construct a B-spline basis matrix and associated penalty matrix for modeling a time variable,
#' including QR decomposition for numerical stability and
#' eigen-decomposition for penalty regularization.
#'
#' @param n_basis Integer. Number of basis functions to use (number of columns in the resulting model matrix).
#' @param time Numeric vector. Values of the time covariate over which the basis is to be constructed.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{model_matrix}{A matrix of dimension length(time) by n_basis. The processed B-spline basis matrix.}
#'     \item{nD}{A square matrix of dimension n_basis by n_basis. The diagonal penalty (or scaling) matrix for regularization/penalization of coefficients.}
#'   }
#'
#' @details
#' Uses degree-3 B-splines with minimum degrees of freedom (df = max(4, n_basis)) including intercept. Applies QR decomposition and transformation for numerical stability. Constructs a second-order difference penalty matrix, computes its eigendecomposition, and separates fixed/random (null space/range space) spline components. The model_matrix output is standardized and aligned in row order with the input time. The penalty matrix nD can be used as a diagonal prior precision or scaling matrix for coefficient penalization.
#'
#' @importFrom splines bs
#' @importFrom stats var
#' @importFrom utils tail
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#'   # Example: Spline basis for irregular time process
#'   time <- sort(sample(1:30, 25, replace = TRUE))
#'   n_basis <- 6
#'   basis <- gen_basis_matrix(n_basis, time)
#' }
#' @export
gen_basis_matrix = function(n_basis, time) {

  order = 1
  B = splines::bs(x = time, df = max(c(4, n_basis)), degree = 3, intercept = TRUE)

  # QR decomposition of B
  qrB = qr(B)
  R = qr.R(qrB)
  Rinv = solve(R)
  Q = qr.Q(qrB)

  # Second-order penalty matrix
  D = diff(diag(max(c(n_basis, 4))), differences = order)
  S = crossprod(D)

  P = t(Rinv) %*% S %*% Rinv

  # Eigen decomposition
  eP = eigen(P)
  evals = eP$values
  evecs = eP$vectors

  # Threshold for numerical zero
  tol = 1e-10
  idx_null = which(evals < tol)
  idx_pos  = which(evals >= tol)

  # Fixed and random components
  U0 = evecs[, idx_null, drop = FALSE]
  Up = evecs[, idx_pos, drop = FALSE]
  lam_pos = evals[idx_pos]

  X = Q %*% U0
  Z = Q %*% Up

  var_z = Z %>% apply(MARGIN = 2, FUN = var)

  model_matrix = cbind(1, scale(Z[, rev(1:ncol(Z))], center = F, scale = T))
  model_matrix = model_matrix[, 1:n_basis]

  rownames(model_matrix) = as.character(time)

  out = list(model_matrix = model_matrix, nD = diag(c(0, rev(evals[idx_pos])/var_z)))

  return(out)
}
