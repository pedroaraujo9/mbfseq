#' Create and Initialize MCMC Sample List for the Cluster Sequence Model
#'
#' Internal: Initialize and return all arrays and lists of parameter values to be
#' sampled or tracked during MCMC fitting, including chains for coefficients, cluster assignments, probabilities, and log-posterior values.
#'
#' @param model_data A model_data object as produced by `create_model_data`,
#' containing all data and model information.
#' @param iters Integer. Number of MCMC iterations for the sample chain (default: 1000).
#' @param init_list (Optional) A list of initial values (for e.g. alpha, w, pw, etc.).
#' If not supplied, values are initialized randomly.
#' @param seed (Optional) Integer random seed for reproducibility.
#'
#' @return A named list containing arrays and matrices to store MCMC samples,
#' including alpha, lambda, w, pw, w_post_prob, id_marg_logp, w_logpost, prob,
#' and (if needed) z, z_post_prob, mu, sigma.
#'
#' @details
#' If init_list is supplied, initial parameter values are set from it; otherwise
#' random initialization is used for MCMC. Dirichlet-distributed cluster weights
#' are sampled with extraDistr.
#' Secondary cluster allocations (w), probabilities, and all arrays are sized according
#' to model data. When model_data$z is absent (i.e., primary clustering is sampled),
#' containers for its samples and cluster parameters (mu, sigma) are initialized.
#' Cluster probabilities are initialized via softmax transformation from mclust.
#' All arrays are preallocated for efficiency and reproducibility.
#'
#' @importFrom stats rnorm var
#' @importFrom extraDistr rdirichlet
#' @importFrom mclust softmax
#' @keywords internal
create_sample_list = function(model_data,
                              iters = 1000,
                              init_list = NULL,
                              seed = NULL) {

  G = model_data$G
  M = model_data$M
  n_basis = model_data$n_basis
  n_id = model_data$n_id
  n_vars = model_data$n_vars
  n = model_data$n
  id_unique = model_data$id_unique
  n_time = model_data$n_time

  if(!is.null(seed)) set.seed(seed)

  alpha = gen_sample_array(
    iters = iters,
    dimension = c(n_basis * M, G),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$alpha
  )

  alpha[,,G] = 0

  if(is.null(model_data$w)) {

    w = gen_sample_array(
      iters = iters,
      dimension = n_id,
      sampler = function(x) sample(1:M, size = n_id, replace = T),
      init = init_list$w
    )

    pw = gen_sample_array(
      iters = iters,
      dimension = c(M),
      sampler = function(x) {
        if(M > 1) {
          extraDistr::rdirichlet(n=1, alpha = rep(10, M)) %>% as.numeric() %>% return()
        }else{
          return(1)
        }
      },

      init = init_list$pw
    )

    w_post_prob = gen_sample_array(
      iters = iters,
      dimension = c(n_id, M),
      init = matrix(pw[1,], nrow = n_id, ncol = M, byrow = T)
    )

  }else{
    w = matrix(model_data$w, nrow = iters, ncol = n_id, byrow = T)
    pw = NULL
    w_post_prob = NULL
  }

  colnames(w) = id_unique

  #### starge probs ####
  X_unique = model_data$X_unique
  stage_prob = gen_sample_array(
    iters = iters,
    dimension = c(n_time*M, G),
    init = mclust::softmax(X_unique %*% alpha[1,, ])
  )

  #### estimating z ####
  if(is.null(model_data$z)){

    z = gen_sample_array(
      iters = iters,
      dimension = n,
      sampler = function(x) sample(1:G, size = n, replace = T),
      init = init_list$z
    )

    z_post_prob = gen_sample_array(
      iters = iters,
      dimension = c(n, G),
      init = 1/G
    )

    mu = gen_sample_array(
      iters = iters,
      dimension = c(G, n_vars),
      sampler = function(x) rnorm(x, sd = 0.1),
      init = init_list$mu
    )

    sigma = gen_sample_array(
      iters = iters,
      dimension = G,
      sampler = function(x) exp(rnorm(x, sd = 0.01)),
      init = init_list$sigma
    )

  }else{
    z = matrix(model_data$z, ncol = n, nrow = iters, byrow = T)
    mu = NULL
    sigma = NULL
    z_post_prob = NULL
  }

  #### out ####
  sample = list(
    alpha = alpha,
    w = w,
    pw = pw,
    w_post_prob = w_post_prob,
    stage_prob = stage_prob,
    z = z,
    z_post_prob = z_post_prob,
    mu = mu,
    sigma = sigma
  )

  return(sample)

}
