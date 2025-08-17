devtools::load_all()
devtools::document()

model_data = create_model_data(
  time = sim_data$data$time,
  id = sim_data$data$id,
  z = sim_data$data$true_z,
  x = NULL,
  w = NULL,
  G = 2,
  M = 2,
  n_basis = 10,
  intercept = FALSE
)

G = 3
M = c(2, 3)
z = sim_data$data$true_z
w = NULL
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

fits = fit_mbfseq(
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

fit = fits

fit$models$`G=3, M=3`$model_info$

fits$metrics
fits$models$
fit = fits$fit$`G=3, M=2`
f %>% posterior_summary()

model = fits$models$`G=3, M=3`

post_summ = fits %>% posterior_summary(M = 3)
post_summ$`G=3, M=3` %>% plot_probability()

fits %>% plot_seq_cluster(M = 2, G = 3)
