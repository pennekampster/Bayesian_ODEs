# Library -----------------------------------------------------------------
library(tidyverse)
library(here)
library(cmdstanr)
## make sure to:
# cmdstanr::install_cmdstan()


# Load and prepare data ---------------------------------------------------------------

D <- read_csv(here("data/FRT_Dataset.csv"))

## returns a list with data structures as declared in "FRT.stan"
getData <- function(D, temp = NA) {
  if(!is.na(temp)) {
    D <- dplyr::filter(D, Temperature == temp)
  }
  data <- list(
    N_series = nrow(D),
    rep_temp = as.integer(as.factor(D$Temperature)),
    time_end = as.matrix(D["Incubation_time"]), # expects array[N_series, 1]
    x_start = cbind(round(D$Prey_start_density), round(D$Predator_start_density)), # expects array[N_series] vector[2]
    x_start_vector = cbind(D$Prey_start_density, D$Predator_start_density), # expects array[N_series] vector[2]
    x_end = cbind(round(D$Prey_end_density), round(D$Predator_end_density))
  )
  return(data)
}


# Fit stan model ----------------------------------------------------------

model <- cmdstan_model("FRT_hierarchical.stan")
n_chains <- 3
# if(!dir.exists("Draws")) dir.create("Draws")
fit <- model$sample(data = getData(D, temp = NA),
                       init = replicate(n_chains, list(b_log = -4.106669, c_log = -5.600194, h_log = -3.531813, K_log = 7.901859, q = 0.867401, r_log = -0.994875), simplify = F), # sigma = 0.222916
                       iter_warmup = 300, iter_sampling = 500, chains = n_chains, parallel_chains = n_chains, output_dir = "Draws", output_basename = "fit_25", seed = 1)

# fit$save_output_files(dir = "Draws", basename = "fit_25")
# fit <- cmdstanr::read_cmdstan_csv(fit_25$output_files())

# Explore fit ----------------------------------------------------------

fit$summary()

includepars <- c("b_log", "c_log","h_log", "K_log", "K_log_temp", "q", "r_log")
draws <- fit$draws()
bayesplot::mcmc_trace(draws, pars = includepars)
bayesplot::mcmc_areas(draws, area_method = "equal height", pars = includepars)
bayesplot::mcmc_pairs(draws, pars = includepars)

