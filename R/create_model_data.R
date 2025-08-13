#' Create Data Structures for Model Building with Dual Cluster Assignments
#'
#' Construct and return a list containing all necessary data structures required for model building based on time, group assignments (z), an optional second clustering (w), and optionally provided covariates and basis functions.
#'
#' @param time Numeric vector indicating the time for each observation.
#' @param id Vector of subject or group identifiers (length = number of observations).
#' @param x (Optional) Matrix or data.frame of covariates (columns = covariates, rows = observations).
#' @param z (Optional) Vector specifying cluster membership for each subject or observation (primary clustering).
#' @param w (Optional) Vector specifying cluster membership for each subject or observation in a second clustering (parallel to z).
#' @param G (Optional) Integer. Number of clusters in z. Required if z is not supplied.
#' @param M (Optional) Integer. Number of clusters in w. Required if w is not supplied.
#' @param n_basis Integer. Number of basis functions for representing time (default is 15).
#' @param intercept Logical. Whether to force an intercept column in cluster matrices (default is FALSE).
#'
#' @return A list of class "model_data" containing prepared basis expansions, cluster encodings, data frames, matrices, and auxiliary structures for downstream model processing.
#'
#' @details
#' Both z and w are interpreted as cluster assignment vectors. If one is not given, you must specify the respective cluster count (G, M). The time variable is expanded with user-defined basis functions (see gen_basis_matrix). Computes auxiliary indices and structures to support penalization and clustering logic for two clusterings. Uses TraMineR for sequence distance calculations among clusters in z. All structures are sized and aligned for downstream modeling.
#'
#' @importFrom dplyr mutate arrange
#' @importFrom tibble tibble
#' @importFrom stats as.dist
#' @importFrom TraMineR seqdef seqdist
#' @examples
#' \dontrun{
#'   # Example with two cluster assignments (z and w) per subject
#'   time <- rep(1:3, 10)
#'   id <- rep(1:10, each = 3)
#'   z <- rep(1:2, each = 15)  # primary clustering
#'   w <- sample(1:2, size=30, replace=TRUE)  # secondary clustering
#'   x <- matrix(rnorm(30), ncol = 1)
#'   data <- create_model_data(time = time, id = id, x = x, z = z, w = w)
#'   str(data)
#' }
#' @export
create_model_data = function(time,
                             id,
                             x = NULL,
                             z = NULL,
                             w = NULL,
                             G = NULL,
                             M = NULL,
                             n_basis = 15,
                             intercept = FALSE) {

  if(!is.null(w)) {
    M = length(unique(w))
  }else{

    if(is.null(M)) {
      stop("Argument `M > 0` should be suplied.")
    }
  }

  if(!is.null(z)) {

    G = length(unique(z))

  }else{

    if(is.null(x)) {
      stop("Data `x` should be suplied when `z` is NULL")
    }

    if(is.null(G)){
      stop("Argument `G > 0` should be suplied.")
    }

  }

  if(!is.null(x)) {

    n_vars = ncol(x)
    col_names = colnames(x)

  }else{

    n_vars = NULL
    col_names = NULL

  }

  order = 1
  time_seq = time - min(time) + 1
  time_unique = unique(time_seq)
  n_time = length(time_unique)
  n = length(time_seq)

  id_unique = unique(id)
  n_id = length(id_unique)

  id_time_df = tibble::tibble(
    id = id,
    time_seq = time_seq,
    time = time
  ) |>
    dplyr::mutate(id = factor(id, levels = id_unique))

  basis_funcions = gen_basis_matrix(
    n_basis = n_basis,
    time = time_seq
  )

  B = basis_funcions$model_matrix
  nD = basis_funcions$nD
  B_unique = B[as.character(time_unique), ]
  B_expand = B[, rep(1:n_basis, times = M)]


  if(intercept == TRUE) {
    w_unique[, 1] = 1
  }

  w_unique = diag(M)
  X_unique = kronecker(w_unique, B_unique)


  w_logp_base = expand.grid(
    id = id_time_df$id |> unique(),
    w = 1:M,
    time_seq = id_time_df$time_seq |> unique()
  ) |>
    dplyr::arrange(w, id, time_seq) |>
    dplyr::mutate(id = factor(id, levels = id_unique))

  if(!is.null(x)) {

    z_logp_base = expand.grid(
      id = id_unique,
      time_seq = id_time_df$time_seq |> unique(),
      z = 1:G
    ) |>
      dplyr::arrange(z, id, time_seq) |>
      dplyr::mutate(id = factor(id, levels = id_unique)) |>
      dplyr::left_join(
        data.frame(x, id = id, time_seq = time_seq), by = c("id", "time_seq")
      )
  }else{
    z_logp_base = NULL
  }

  notpen_index = gen_notpen_index(n_basis = n_basis, M = M, order = order)
  basis_index = gen_basis_index(n_basis = n_basis, M = M)

  out = list(
    G = G,
    M = M,
    n_id = n_id,
    n_time = n_time,
    n_basis = n_basis,
    n = n,
    n_vars = n_vars,
    order = order,
    z = z,
    w = w,
    x = x,
    id_time_df = id_time_df,
    time_unique = time_unique,
    id_unique = id_unique,
    col_names = col_names,
    B = B,
    B_unique = B_unique,
    B_expand = B_expand,
    X_unique = X_unique,
    notpen_index = notpen_index,
    basis_index = basis_index,
    w_logp_base = w_logp_base,
    z_logp_base = z_logp_base,
    nD = nD
  )

  class(out) = "model_data"
  return(out)

}
