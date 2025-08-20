# mbfseq: Clustering Sequence Data with Mixtures of Spline Categorical Data

## Overview

The `mbfseq` package provides tools for clustering sequence data using a novel Bayesian approach based on mixtures of spline categorical distributions. This method is particularly suited for sequences observed over time or another continuous covariate, allowing for flexible modelling of the probability of observing different categories along the sequence.

## Key Features

- **Flexible Bayesian clustering** of categorical sequences
- **Spline-based modelling** for smooth probability curves over time/covariates  
- **MCMC algorithms** for robust parameter estimation
- **Model selection** across different numbers of clusters (G) and mixture components (M)
- **Visualization tools** for exploring clustering results
- **Parallel processing** support for faster computation

## Installation

You can install the development version of mbfseq from GitHub:

```r
# Install devtools if you haven't already
if (!require(devtools)) {
  install.packages("devtools")
}

# Install mbfseq
devtools::install_github("pedroaraujo9/mbfseq")
```

## Quick Start

```r
library(mbfseq)

# Load example data
data(sim_data)

# Fit the model with different cluster numbers
result <- fit_mbfseq(
  G = c(2, 3),           # Number of clusters to try
  M = c(3, 4),           # Number of mixture components to try
  id = sim_data$id,      # Subject identifiers
  time = sim_data$time,  # Time points
  iters = 1000,          # MCMC iterations
  n_cores = 2            # Parallel processing
)

# View model metrics
result$metrics

# Plot clustering results
plot_seq_cluster(fit = result, G = 2, M = 3)

# Plot probability curves
plot_probability(fit = result, G = 2, M = 3)

# Get posterior summaries
posterior_summary(result, G = 2, M = 3)
```

## Main Functions

### Model Fitting
- `fit_mbfseq()`: Main function for fitting the Bayesian mixture model

### Visualization
- `plot_seq_cluster()`: Visualize sequence clustering as heatmap
- `plot_probability()`: Plot probability curves over time



## Example: Simulated Data

```r
# Generate simulated data
set.seed(123)
n_subjects <- 50
n_timepoints <- 10

# Create example sequence data
sim_data <- data.frame(
  id = rep(1:n_subjects, each = n_timepoints),
  time = rep(1:n_timepoints, times = n_subjects),
  category = sample(1:4, n_subjects * n_timepoints, replace = TRUE)
)

# Fit model
fit <- fit_mbfseq(
  G = 2:3,
  M = 3:4,
  id = sim_data$id,
  time = sim_data$time,
  z = sim_data$category,  # Optional: provide sequence data directly
  iters = 500,
  burn_in = 250,
  thin = 2
)

# Examine results
print(fit$metrics)
plot_seq_cluster(fit = fit, G = 2, M = 3)
```

## Advanced Usage

### Custom Configuration

```r
# Advanced configuration
config <- list(
  bounds = c(0.01, 10),      # Lambda search bounds
  n_start = 50,              # Number of random initializations
  n_start_iters = 30,        # Iterations for initialization search
  epsilon_w = 1,             # Prior parameters
  beta_sd = sqrt(10),
  mu_sd = sqrt(10),
  sigma_a = 1,
  sigma_b = 1
)

fit <- fit_mbfseq(
  G = 2:4,
  M = 3:5,
  id = sim_data$id,
  time = sim_data$time,
  config = config,
  verbose = TRUE
)
```

## Dependencies

The package depends on several R packages:
- `Rcpp` and `RcppArmadillo` for efficient computation
- `ggplot2` and `viridis` for visualization
- `dplyr` and `tidyr` for data manipulation
- `future` and `future.apply` for parallel processing
- `TraMineR` for sequence analysis utilities


## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests on [GitHub](https://github.com/pedroaraujo9/mbfseq).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


