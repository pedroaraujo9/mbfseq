#' Run a Single MCMC Sampling Chain for the Model
#'
#' This function runs a single Markov Chain Monte Carlo (MCMC) sampling chain
#' for the specified model using Gibbs updates for latent variables and parameters.
#' It supports optional verbose output and configurable priors, iterations, burn-in, and thinning.
#'
#' @param model_data A list containing the model data and settings required for sampling.
#' @param lambda Numeric scalar controlling regularization strength (default 1).
#' @param iters Integer number of MCMC iterations (default 1000).
#' @param burn_in Number of iterations to discard as burn-in (default half of \code{iters}).
#' @param thin Integer thinning interval for stored samples (default 10).
#' @param init_list Optional list of initial parameter values (default \code{NULL}).
#' @param priors A list of prior hyperparameters:
#' \describe{
#'   \item{epsilon_w}{Prior parameter for updating \code{w} (default 1).}
#'   \item{beta_sd}{Prior standard deviation for \code{alpha} coefficients (default \code{sqrt(10)}).}
#'   \item{mu_sd}{Prior standard deviation for \code{mu} parameters (default \code{sqrt(10)}).}
#'   \item{sigma_a}{Shape parameter for \code{sigma} prior (default 1).}
#'   \item{sigma_b}{Rate parameter for \code{sigma} prior (default 1).}
#' }
#' @param verbose Logical indicating whether to print progress updates (default \code{FALSE}).
#' @param seed Optional integer seed for reproducibility (default \code{NULL}).
#'
#' @return A list of class \code{"mbfseq"} containing:
#' \itemize{
#'   \item \code{z_class}: Estimated latent class assignments for \code{z}.
#'   \item \code{w_class}: Estimated latent class assignments for \code{w}.
#'   \item \code{sample_list}: A list of sampled parameter arrays.
#'   \item \code{seed}: The seed used when fitting the model.
#' }
#'
#' @details
#' The function performs Gibbs updates for latent parameters \code{alpha}, \code{w}, \code{pw}, \code{z}, \code{mu}, and \code{sigma}
#' iteratively over \code{iters} MCMC steps. It supports burn-in and thinning to reduce autocorrelation in stored samples.
#' The updates depend on the provided priors and initial values.
#' Progress updates are printed every 10% of iterations if \code{verbose = TRUE}.
#'
#' @importFrom purrr map
#' @importFrom utils flush.console
#'
#' @examples
#' \dontrun{
#' # Assuming 'model_data' is prepared and required dependencies are loaded
#' fit <- single_run(model_data, iters = 2000, verbose = TRUE, seed = 123)
#' summary(fit$z_class)
#' }
single_run = function(model_data,
                      lambda = 1,
                      iters = 1000,
                      burn_in = iters/2,
                      thin = 10,
                      init_list = NULL,
                      priors = list(
                        epsilon_w = 1,
                        beta_sd = sqrt(10),
                        mu_sd = sqrt(10),
                        sigma_a = 1,
                        sigma_b = 1
                      ),
                      verbose = FALSE,
                      seed = NULL) {

  if(is.null(seed)) seed = sample(1:10000, size = 1)
  set.seed(seed)

  M = model_data$M
  G = model_data$G
  intercept = FALSE

  sample_list = create_sample_list(
    iters = iters,
    model_data = model_data,
    init_list = init_list,
    seed = NULL
  )

  logpost_init = eval_logpost(
    z = sample_list$z[1,],
    w = sample_list$w[1, ],
    alpha = sample_list$alpha[1,,],
    mu = sample_list$mu[1,,],
    sigma = sample_list$sigma[1,],
    lambda = lambda,
    priors = priors,
    model_data = model_data
  )

  logpost = gen_sample_array(
    iters = iters,
    dimension = 3,
    init = logpost_init
  )

  colnames(logpost) = names(logpost_init)

  print_h = floor(iters/10)

  for(i in 2:iters) {

    if(verbose) {
      if (i %% print_h == 0) {
        cat(sprintf("M = %d - Iteration %d / %d\n", M, i, iters))
        utils::flush.console()
      }
    }

    #### update alpha, w, pw ####
    sample_list$alpha[i,,] = update_alpha(
      alpha = sample_list$alpha[i-1,,],
      lambda = lambda,
      z = sample_list$z[i-1, ],
      w = sample_list$w[i-1, ],
      model_data = model_data,
      fixed_sd = priors$beta_sd,
      intercept = intercept
    )

    sample_list$stage_prob[i,,] = comp_prob(
      alpha = sample_list$alpha[i,,],
      model_data = model_data
    )

    if(is.null(model_data$w)) {

      sample_list$pw[i, ] = update_pw(
        w = sample_list$w[i-1,],
        epsilon = priors$epsilon_w,
        model_data = model_data
      )

      new_w = update_w(
        alpha = sample_list$alpha[i,,],
        pw = sample_list$pw[i, ],
        z = sample_list$z[i-1, ],
        model_data = model_data
      )

      sample_list$w[i, ] = new_w$w
      sample_list$w_post_prob[i,,] = new_w$w_post_prob

    }

    #### update z, mu, sigma ####

    if(is.null(model_data$z)) {
      sample_list$mu[i,,] = update_mu(
        z = sample_list$z[i-1, ],
        sigma = sample_list$sigma[i-1, ],
        model_data = model_data,
        mu_fixed_sd = priors$mu_sd
      )

      sample_list$sigma[i, ] = update_sigma(
        mu = sample_list$mu[i,,],
        z = sample_list$z[i-1, ],
        sigma_a = priors$sigma_a,
        sigma_b = priors$sigma_b,
        model_data = model_data
      )

      z_up = update_z(
        mu = sample_list$mu[i,,],
        sigma = sample_list$sigma[i, ],
        w = sample_list$w[i,],
        alpha = sample_list$alpha[i,,],
        model_data = model_data
      )

      sample_list$z[i, ] = z_up$z
      sample_list$z_post_prob[i,,] = z_up$z_post_prob
    }

    logpost[i, ] = eval_logpost(
      z = sample_list$z[i,],
      w = sample_list$w[i, ],
      alpha = sample_list$alpha[i,,],
      mu = sample_list$mu[i,,],
      sigma = sample_list$sigma[i,],
      lambda = lambda,
      priors = priors,
      model_data = model_data
    )
  }

  if(verbose) {
    cat("\n")
  }

  sample_list = purrr::map(sample_list, filter_array, burn_in = burn_in, thin = thin)
  logpost = filter_array(logpost, burn_in = burn_in, thin = thin)

  out = list(
    z_class = sample_list$z |> comp_class(),
    w_class = sample_list$w |> comp_class(),
    sample_list = sample_list,
    logpost = logpost,
    seed = seed
  )

  class(out) = "mbfseq"
  return(out)
}
