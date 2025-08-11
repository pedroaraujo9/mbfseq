devtools::load_all()

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

