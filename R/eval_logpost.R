#' Evaluate Full Log-Posterior for Sequence Clustering Model
#'
#' Computes the complete log-posterior density for a sequence clustering model, combining:
#' \itemize{
#'   \item Sequence model likelihood (via \code{eval_seq_logpost})
#'   \item Observation model likelihood (Gaussian)
#'   \item Prior distributions for all parameters
#' }
#'
#' @param z Integer vector of cluster assignments (latent sequence memberships)
#' @param w Numeric vector of cluster-specific latent variables
#' @param alpha Numeric matrix of category-specific coefficients
#' @param mu Numeric matrix of cluster means (G x p)
#' @param sigma Numeric vector of cluster standard deviations
#' @param lambda Numeric penalty parameter for regularization
#' @param priors List containing prior parameters:
#' \describe{
#'   \item{epsilon_w}{Precision for w prior (default: 1)}
#'   \item{beta_sd}{SD for non-penalized coefficients (default: sqrt(10))}
#'   \item{mu_sd}{SD for cluster means prior (default: sqrt(10))}
#'   \item{sigma_a}{Shape parameter for sigma prior (default: 1)}
#'   \item{sigma_b}{Rate parameter for sigma prior (default: 1)}
#' }
#' @param model_data List containing model structure and data components
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{logpost}{Complete log-posterior including all terms}
#'   \item{penal_logpost}{Log-posterior including only penalized terms}
#'   \item{logp_seq_model}{Log-posterior for just the sequence model}
#' }
#'
#' @details
#' The function computes:
#' \enumerate{
#'   \item Sequence model log-posterior using \code{eval_seq_logpost}
#'   \item If \code{model_data$z} is NULL, adds:
#'   \itemize{
#'     \item Gaussian likelihood for observations x
#'     \item Gaussian prior for cluster means mu
#'     \item Inverse-gamma prior for cluster variances sigma^2
#'   }
#' }
#'
#' The inverse-gamma prior is parameterized such that:
#' sigma^2 ~ InvGamma(shape = sigma_a, rate = sigma_b)
#'
#' @importFrom stats dnorm
#' @importFrom extraDistr dinvgamma
#'
#' @examples
#' \dontrun{
#' # Example setup (not run)
#' data <- list(
#'   x = matrix(rnorm(100*3), ncol=3),
#'   G = 3,
#'   n = 100,
#'   z = NULL,
#'   # ... other required model_data components ...
#' )
#'
#' priors <- list(
#'   epsilon_w = 1,
#'   beta_sd = sqrt(10),
#'   mu_sd = sqrt(10),
#'   sigma_a = 1,
#'   sigma_b = 1
#' )
#'
#' result <- eval_logpost(
#'   z = sample(1:3, 100, replace=TRUE),
#'   w = rnorm(3),
#'   alpha = matrix(rnorm(6), nrow=2, ncol=3),
#'   mu = matrix(rnorm(9), nrow=3, ncol=3),
#'   sigma = runif(3),
#'   lambda = 1.0,
#'   priors = priors,
#'   model_data = data
#' )
#' }
#'
#' @seealso \code{\link{eval_seq_logpost}} for the sequence model component
#' @export
eval_logpost = function(z,
                        w,
                        alpha,
                        mu,
                        sigma,
                        lambda,
                        priors = list(
                          epsilon_w = 1,
                          beta_sd = sqrt(10),
                          mu_sd = sqrt(10),
                          sigma_a = 1,
                          sigma_b = 1
                        ),
                        model_data) {
  x = model_data$x
  G = model_data$G
  n = model_data$n

  logp_seq_model = eval_seq_logpost(
    alpha = alpha,
    w = w,
    z = z,
    lambda = lambda,
    model_data = model_data,
    fixed_sd = priors$beta_sd
  )

  if(is.null(model_data$z)) {

    mu_expand = mu[z, ]
    sigma_expand =  matrix(sigma[z], nrow = n, ncol = G)

    logp_x = sum(dnorm(x, mean = mu_expand, sd = sigma_expand, log = TRUE))
    logp_mu = sum(dnorm(mu, mean = 0, sd = priors$mu_sd, log = TRUE))

    logp_sigma = sum(
      extraDistr::dinvgamma(sigma^2, alpha = priors$sigma_a, beta = priors$sigma_b, log = TRUE)
    )

  }else{

    logp_x = 0
    logp_mu = 0
    logp_sigma = 0

  }


  logpost = logp_seq_model + logp_x + logp_mu + logp_sigma
  penal_logpost = logp_seq_model + logp_x

  out = c(
    logpost = logpost,
    penal_logpost = penal_logpost,
    logp_seq_model = logp_seq_model
  )

  return(out)

}
