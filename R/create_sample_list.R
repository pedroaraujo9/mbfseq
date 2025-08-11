#' Create and Initialize MCMC Sample List for the Cluster Sequence Model
#'
#' Internal: Initialize and return all arrays and lists of parameter values to be sampled or tracked during MCMC fitting, including chains for coefficients, cluster assignments, probabilities, and log-posterior values.
#'
#' @param model_data A model_data object as produced by create_model_data, containing the full model structure and cluster counts.
#' @param iters Integer. Number of MCMC iterations for the sample chain (default: 1000).
#' @param epsilon Numeric. Small value used for Dirichlet prior on cluster weights (default: 1).
#' @param lambda (Optional) Value(s) for the regularization parameter(s); can be a scalar or matrix as used in the model.
#' @param init_list (Optional) A list of initial values (for e.g. alpha, w, pw, etc.). If not supplied, values are initialized randomly.
#' @param seed (Optional) Integer random seed for reproducibility.
#'
#' @return A named list containing arrays and matrices to store MCMC samples, including alpha, lambda, w, pw, w_post_prob, id_marg_logp, w_logpost, prob, and (if needed) z, z_post_prob, mu, sigma.
#'
#' @details
#' If init_list is supplied, initial parameter values are set from it; otherwise random initialization is used for MCMC. Dirichlet-distributed cluster weights are sampled with extraDistr. Secondary cluster allocations (w), probabilities, and all arrays are sized according to model data. When model_data$z is absent (i.e., primary clustering is sampled), containers for its samples and cluster parameters (mu, sigma) are initialized. Cluster probabilities are initialized via softmax transformation from mclust. All arrays are preallocated for efficiency and reproducibility.
#'
#' @importFrom stats rnorm var
#' @importFrom extraDistr rdirichlet
#' @importFrom mclust softmax
#' @keywords internal
create_sample_list = function(model_data,
                              iters = 1000,
                              epsilon = 1,
                              lambda = NULL,
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

  if(!is.null(init_list)) {
    alpha_init = init_list$alpha
    w_init = init_list$w
    pw_init = init_list$pw
    z_init = init_list$z_init
    mu_init = init_list$mu_init
    sigma_init = init_list$sigma_init

  }else{

    alpha_init = matrix(
      stats::rnorm(n_basis * M * G, sd = 0.01),
      nrow = n_basis * M,
      ncol = G
    )

    w_init = sample(1:M, size = n_id, replace = T)

    if(M > 1) {

      pw_init = extraDistr::rdirichlet(n=1, alpha = rep(epsilon + 100, M)) |> as.numeric()

    }else{

      pw_init = 1

    }

    if(is.null(model_data$z)) {

      z_init = sample(1:G, size = n, replace = T)

      mu_init = matrix(
        stats::rnorm(G * n_vars, sd = 0.01), nrow = G, ncol = n_vars
      )

      sigma_init = exp(stats::rnorm(G, sd = 0.01))

    }
  }

  #### alpha and lambda ####
  alpha = array(0, dim = c(iters, n_basis * M, G))
  alpha[1, ,] = alpha_init
  alpha[1,, G] = 0

  lambda = array(lambda, dim = c(iters, M, G-1))

  #### w ####
  w = array(NA, dim = c(iters, n_id))
  colnames(w) = id_unique

  if(is.null(model_data$w)) {

    w_post_prob = array(0, dim = c(iters, n_id, M))
    id_marg_logp = array(0, dim = c(iters, n_id))
    pw = array(0, dim = c(iters, M))

    w[1, ] = w_init
    pw[1, ] = pw_init
    w_post_prob[1, , ] = 1/M
    id_marg_logp[1, ] = 0

  }else{

    w = matrix(model_data$w, nrow = iters, ncol = n_id, byrow = T)
    colnames(w) = id_unique
    pw = NULL
    w_post_prob = NULL
    id_marg_logp  = NULL

  }

  #### z ####
  if(is.null(model_data$z)) {
    mu = array(NA, dim = c(iters, G, n_vars))
    sigma = array(NA, dim = c(iters, G))

    z = array(NA, dim = c(iters, n))
    z_post_prob = array(NA, dim = c(iters, n, G))

    mu[1,,] = mu_init
    sigma[1, ] = sigma_init

    z[1, ] = z_init
    z_post_prob[1,,] = 1/G

  }else{
    mu = NULL
    sigma = NULL
    z = NULL
    z_post_prob = NULL
  }

  #### probs ####
  X_unique = model_data$X_unique
  prob = array(NA, dim = c(iters, n_time*M, G))
  prob[1, , ] = mclust::softmax(X_unique %*% alpha[1,, ])

  #### logpost ####
  w_logpost = numeric(length = iters)
  z_logpost = numeric(length = iters)

  #### out ####
  sample = list(
    alpha = alpha,
    lambda = lambda,
    w = w,
    pw = pw,
    w_post_prob = w_post_prob,
    id_marg_logp  = id_marg_logp,
    w_logpost = w_logpost,
    prob = prob,
    z = z,
    z_post_prob = z_post_prob,
    mu = mu,
    sigma = sigma,
    z_logpost = z_logpost
  )

  return(sample)

}
