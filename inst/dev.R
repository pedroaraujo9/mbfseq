devtools::load_all()
devtools::document()


G = 3
M = c(2, 3)
z = sim_data$data$true_z %>% factor(labels = c("A", "B", "C"))
w = rep(1, length(sim_data$true_w))
x = NULL
id = sim_data$data$id
time = sim_data$data$time
iters = 100
burn_in = iters/2
thin = 2
lambda = NULL
n_basis = 10
init_list = NULL
n_cores = 1
config = list(
  bounds = c(0.01, 5),
  lambda_start = 1,
  n_points = 20,
  n_start = 20,
  n_start_iters = 10,
  n_start_cores = 1,
  epsilon_w = 1,
  beta_sd = sqrt(10),
  mu_sd = sqrt(10),
  sigma_a = 1,
  sigma_b = 1
)
verbose = TRUE
seed = NULL

fit = fit_mbfseq(
  G = G,
  M = M,
  z = z,
  w = w,
  x = x,
  id = id,
  time = time,
  iters = iters,
  burn_in = burn_in,
  thin = thin,
  lambda = lambda,
  n_basis = n_basis,
  init_list = init_list,
  n_cores = n_cores,
  config = config,
  verbose = verbose,
  seed = seed
)

fit$models$`G=3, M=2`$sample_list$w
fit$models$`G=3, M=3`$model_info$

fit$
fits$models$
fit = fits$fit$`G=3, M=2`
f %>% posterior_summary()

model = fits$models$`G=3, M=3`

post_summ = fits %>% posterior_summary(M = 3)
post_summ$`G=3, M=3` %>% plot_probability()

fit$model_data$z_levels


devtools::load_all()
fit %>% plot_seq_cluster(M = 2, G = 3)
fit %>% plot_seq_cluster(M = 3, G = 3)


