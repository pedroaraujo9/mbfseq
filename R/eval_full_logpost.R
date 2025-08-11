eval_full_logpost = function(z,
                             w,
                             alpha,
                             mu,
                             sigma,
                             lambda,
                             spline_fixed_sd,
                             mu_fixed_sd,
                             sigma_a,
                             sigma_b,
                             model_data) {
  x = model_data$x
  G = model_data$G
  n = model_data$n

  mu_expand = mu[z, ]
  sigma_expand = sigma[z] |> matrix(nrow = n, ncol = G)

  logp_xz = sum(dnorm(x, mean = mu_expand, sd = sigma_expand, log = TRUE))
  logp_mu = sum(dnorm(mu, mean = 0, sd = mu_fixed_sd, log = TRUE))
  logp_sigma = sum(
    extraDistr::dinvgamma(sigma^2, alpha = sigma_a, beta = sigma_b, log = TRUE)
  )

  logp_zw = comp_logpost(
    alpha = alpha,
    w = w,
    z = z,
    lambda = lambda,
    model_data = model_data,
    fixed_sd = spline_fixed_sd
  )

  logpost = logp_xz + logp_zw + logp_mu + logp_sigma
  logplike = logp_xz + logp_zw

  out = c(
    logpost = logpost,
    logplike = logplike,
    logp_xz = logp_xz,
    logp_zw = logp_zw,
    logp_mu = logp_mu,
    logp_sigma = logp_sigma
  )

  return(out)

}
