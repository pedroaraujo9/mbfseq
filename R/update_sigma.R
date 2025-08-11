update_sigma = function(mu, z, sigma_a = 1, sigma_b = 1, model_data) {
  x = model_data$x
  G = model_data$G
  n_vars = model_data$n_vars

  sigma2 = numeric(G)

  for (g in 1:G) {
    idx = which(z == g)
    n_k = length(idx)

    x_k = x[idx, , drop = FALSE]
    mu_k = matrix(mu[g, ], nrow = n_k, ncol = n_vars, byrow = TRUE)

    sse = sum((x_k - mu_k)^2)

    shape = sigma_a + (n_k * n_vars) / 2
    rate = sigma_b + sse / 2

    sigma2[g] = 1 / rgamma(1, shape = shape, rate = rate)
  }

  sigma = sqrt(sigma2)
  return(sigma)
}
