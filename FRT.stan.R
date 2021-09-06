# Library -----------------------------------------------------------------
library(tidyverse)
library(here)
library(cmdstanr)
library(here)

library(posterior)
library(bayesplot)
color_scheme_set("brightblue")


## make sure to:
# remotes::install_github("stan-dev/cmdstanr")
# cmdstanr::install_cmdstan()


# Load and prepare data ---------------------------------------------------------------

D <- read_csv(here("data/FRT_Dataset.csv"))

## returns a list with data structures as declared in "FRT.stan"
getData <- function(D, temp = 25) {
  D <- dplyr::filter(D, Temperature == temp & Treatment == "predator")
  data <- list(
    N_series = nrow(D),
    time_end = as.matrix(D["Incubation_time"]), # expects array[N_series, 1]
    x_start = cbind(round(D$Prey_start_density), round(D$Predator_start_density)), # expects array[N_series] vector[2]
    x_start_vector = cbind(D$Prey_start_density, D$Predator_start_density), # expects array[N_series] vector[2]
    x_end = cbind(round(D$Prey_end_density), round(D$Predator_end_density))
  )
  return(data)
}


# Fit stan model ----------------------------------------------------------

model <- cmdstan_model("FRT.stan")
n_chains <- 3
# if(!dir.exists("Draws")) dir.create("Draws")
fit_25 <- model$sample(data = getData(D, temp = 25),
                       init = replicate(n_chains, list(b_log = -4.106669, c_log = -5.600194, h_log = -3.531813, K_log = 7.901859, q = 0.867401, r_log = -0.994875), simplify = F), # sigma = 0.222916
                       iter_warmup = 300, iter_sampling = 300, chains = n_chains, parallel_chains = n_chains, output_dir = "Draws", output_basename = "fit_25", seed = 1)

# fit_25$save_output_files(dir = "Draws", basename = "fit_25")
# fit_25 <- cmdstanr::read_cmdstan_csv(fit_25$output_files())

# Explore fit ----------------------------------------------------------

fit_25$summary()
includepars <- c("b_log", "c_log","h_log", "K_log", "q", "r_log")
draws_25 <- fit_25$draws()
bayesplot::mcmc_trace(draws_25, pars = includepars)
bayesplot::mcmc_areas(draws_25, area_method = "equal height", pars = includepars)
bayesplot::mcmc_pairs(draws_25, pars = includepars)

