update_mu = function(z, sigma, model_data, mu_fixed_sd = 1) {

  G = model_data$G
  x = model_data$x
  n_vars = model_data$n_vars

  mu = matrix(NA, nrow = G, ncol = n_vars)
  ng = z %>% factor(levels = 1:G) %>% table()

  for (g in 1:G) {
    var_g = 1/((1/(mu_fixed_sd^2)) + ng[g]/(sigma[g]^2))
    mu_g = var_g * (ng[g] * colMeans(x[z == g, ])/(sigma[g]^2))
    mu[g, ] = rnorm(n_vars, mean = mu_g, sd = sqrt(var_g))
  }

  return(mu)

}
