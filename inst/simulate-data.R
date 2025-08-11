library(tidyverse)

set.seed(1)
n_id = 30
n_time = 50
id_label = str_c("id", 1:n_id)
time_label = 1:n_time
M = 3
G = 3
wtrue = vector(length = n_id)
wtrue = sample(1:M, size = n_id, replace = T, prob = c(0.4, 0.3, 0.4))
names(wtrue) = id_label

id_effect = rnorm(n = n_id * G, sd = 1) %>% matrix(nrow = n_id, ncol = G)
id_effect[, G] = 0
id_effect

#### generating z ####

# M == 1
m1 = cbind(
  -0.1*(time_label - 10),
  0.5*(time_label - 40),
  0
)

m1 = m1/0.7

p1 = mclust::softmax(m1)

# M == 2
m2 = cbind(
  -2 -scale((time_label - 20)^2),
  scale((time_label - 20)^2),
  0
)

m2 = m2/0.5

p2 = mclust::softmax(m2)

data.frame(time_label, p2) %>%
  gather(group, value, -time_label) %>%
  ggplot(aes(x=time_label, y=value, color=group)) +
  geom_line()

# M == 3
m3 = cbind(
  -2 + scale(sin(2 * pi * time_label / 30)),
  -0 + scale(-time_label^2),
  0
)

m3 = m3

p3 = mclust::softmax(m3)

data.frame(time_label, p3) %>%
  gather(group, value, -time_label) %>%
  ggplot(aes(x=time_label, y=value, color=group)) +
  geom_line()

true_effect = list(m1, m2, m3)

true_prob = lapply(1:M, function(m){
  true_effect[[m]] %>%
    mclust::softmax() %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(time = 1:n_time, w = m) %>%
    gather(cat, prob, -time, -w) %>%
    mutate(cat = str_replace(cat, "V", ""))
}) %>%
  do.call(rbind, .)

true_prob %>%
  ggplot(aes(x=time, y=prob, color=cat)) +
  geom_line() +
  facet_grid(. ~ w)

set.seed(1)
sim_data = lapply(1:n_id, function(i) {

  i_id_effect = cbind(rep(1, n_time)) %*% rbind(id_effect[i, ])
  data.frame(
    id = id_label[i],
    true_z = extraDistr::rcat(n = n_time, prob = mclust::softmax(true_effect[[wtrue[i]]])),
    true_z_id_effect = extraDistr::rcat(n = n_time, prob = mclust::softmax(true_effect[[wtrue[i]]]+i_id_effect)),
    time = time_label,
    true_w = wtrue[i]
  )
}) %>% do.call(rbind, .) %>%
  as_tibble()

true_prob = lapply(1:M, function(m){
  true_effect[[m]] %>%
    softmax() %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(time = 1:n_time, w = m) %>%
    gather(cat, prob, -time, -w) %>%
    mutate(cat = str_replace(cat, "V", ""))
}) %>%
  do.call(rbind, .)

sim_data = list(
  data = sim_data,
  probs = true_prob,
  true_w = wtrue,
  id_effect = id_effect
)

sim_data$data %>%
  mutate(id = factor(id, levels = names(sort(wtrue)))) %>%
  ggplot(aes(x=time, y=id, fill=factor(true_z))) +
  geom_tile(color="white") +
  viridis::scale_fill_viridis(discrete = T)

sim_data$data %>%
  mutate(id = factor(id, levels = names(sort(wtrue)))) %>%
  ggplot(aes(x=time, y=id, fill=factor(true_z_id_effect))) +
  geom_tile(color="white") +
  viridis::scale_fill_viridis(discrete = T)

#### generating x ####
mu = rbind(
  c(1, 2),
  c(-2, -2),
  c(-3, -3)
)

cov_list = list(
  rbind(
    c(3, -2),
    c(-2, 3)
  ),

  rbind(
    c(1, 0.4),
    c(0.4, 1)
  ),

  rbind(
    c(5, 1),
    c(1, 5)
  )
)

mu

cov_list


n = length(sim_data$data$true_z)
z = sim_data$data$true_z
x = matrix(nrow = n, ncol = 2)
sim_data$true_mu = mu

for(i in 1:n) {
  x[i, ] = mvtnorm::rmvnorm(n = 1, mean = mu[z[i], ], sigma = cov_list[[z[i]]])
}

colnames(x) = paste0("var", 1:2)

data.frame(x) %>%
  as_tibble() %>%
  mutate(z = z) %>%
  ggplot(aes(x=var1, y=var2, color=factor(z))) +
  geom_point()

sim_data$x = x



z_random = sample(1:G, size = n, replace = T)
x_random = matrix(nrow = n, ncol = 2)
for(i in 1:n) {
  x_random[i, ] = mvtnorm::rmvnorm(n = 1, mean = mu[z_random[i], ], sigma = diag(2))
}

colnames(x_random) = paste0("var", 1:2)

data.frame(x_random) %>%
  as_tibble() %>%
  mutate(z = z_random) %>%
  ggplot(aes(x=var1, y=var2, color=factor(z))) +
  geom_point()

sim_data$data$z_random = z_random
sim_data$x_random = x_random
#### saving ####
#usethis::use_data(sim_data, internal = TRUE, overwrite = TRUE)
usethis::use_data(sim_data, overwrite = TRUE)





