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

fits$metrics
fits$models$
fit = fits$fit$`G=3, M=2`
f %>% posterior_summary()

model = fits$models$`G=3, M=3`

post_summ = fits %>% posterior_summary(M = 3)
post_summ$`G=3, M=3` %>% plot_probability()

devtools::load_all()
fit %>% plot_seq_cluster(M = 2, G = 3)
fit %>% plot_seq_cluster(M = 3, G = 3)


# Example heatmap data
df <- expand.grid(
  time = 1:10,
  country = c("A", "B", "C", "D")
)
df$category <- sample(c("cat1", "cat2", "cat3"), nrow(df), replace = TRUE)

# Cluster info for each country
clusters <- data.frame(
  country = c("A", "B", "C", "D"),
  cluster = c("Cluster 1", "Cluster 1", "Cluster 2", "Cluster 2")
)
p <- ggplot(df, aes(x = time, y = country, fill = category)) +
  geom_tile(color = "white") +
  scale_fill_brewer(palette = "Set3")

df_labels <- df %>%
  select(country) %>%
  distinct() %>%
  left_join(clusters, by = "country")

# Add text on right side
p + geom_text(
  data = df_labels,
  aes(x = max(df$time) + 0.5, y = country, label = cluster),
  inherit.aes = FALSE,
  hjust = 0
) +
  coord_cartesian(clip = 'off') + # allow text outside plot
  theme(
    plot.margin = margin(5, 50, 5, 5) # extra space on right
  )

