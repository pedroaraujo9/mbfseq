#' Update Cluster Membership Probabilities
#'
#' Internal: Update the latent cluster membership probabilities (pw) using a Dirichlet posterior, based on current cluster allocations.
#'
#' @param w Integer vector of latent cluster memberships (values from 1 to M).
#' @param epsilon Numeric. Value added to Dirichlet concentration parameters (default 1).
#' @param model_data List. Model data structure, must include number of clusters M.
#'
#' @return Numeric vector of updated cluster probabilities (length M).
#'
#' @details
#' This function computes the posterior Dirichlet parameters for the cluster membership probabilities by counting the number of assignments to each cluster (from w), adding epsilon for smoothing, and sampling from the resulting Dirichlet distribution. Probabilities smaller than 1e-300 are thresholded for numerical stability.
#'
#' @importFrom extraDistr rdirichlet
#' @keywords internal
update_pw = function(w, epsilon = 1, model_data) {
  M = model_data$M
  if(M > 1) {
    w = factor(w, levels = 1:M, ordered = T)
    nw = w |> table()
    pw = extraDistr::rdirichlet(1, alpha = as.numeric(nw + epsilon)) |> as.numeric()
    pw[pw < 1e-300] = 1e-300
  }else{
    pw = 1
  }

  return(pw)
}
