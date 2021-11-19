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
  D <- dplyr::filter(D, Temperature == temp)
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
fit_15 <- model$sample(data = getData(D, temp = 15),
                       init = replicate(n_chains, list(b_log = -4.106669, c_log = -5.600194, h_log = -3.531813, K_log = 7.901859, q = 0.867401, r_log = -0.994875), simplify = F), # sigma = 0.222916
                       iter_warmup = 300, iter_sampling = 500, chains = n_chains, parallel_chains = n_chains, output_dir = "Draws", output_basename = "fit_15", seed = 1)

# if(!dir.exists("Draws")) dir.create("Draws")
fit_20 <- model$sample(data = getData(D, temp = 20),
                       init = replicate(n_chains, list(b_log = -15.229665, c_log = -4.554020, h_log = -3.100683, K_log = 7.897008, q = 3.492671, r_log = 0.136438), simplify = F), # sigma = 0.222916
                       iter_warmup = 300, iter_sampling = 500, chains = n_chains, parallel_chains = n_chains, output_dir = "Draws", output_basename = "fit_20", seed = 1)

# if(!dir.exists("Draws")) dir.create("Draws")
fit_25 <- model$sample(data = getData(D, temp = 25),
                       init = replicate(n_chains, list(b_log = 1.117, c_log = -4.60517, h_log = -6.907755, K_log = 8.477, q = -.56, r_log = 0.1897936), simplify = F), # sigma = 0.222916
                       iter_warmup = 300, iter_sampling = 500, chains = n_chains, parallel_chains = n_chains, output_dir = "Draws", output_basename = "fit_25", seed = 1)



#fit_25$save_output_files(dir = "Draws", basename = "fit_25")
#fit_25 <- cmdstanr::read_cmdstan_csv(fs::dir_ls(here("Draws/"), regexp = "fit"))



fit_15 <- as_cmdstan_fit(
  fs::dir_ls("Draws/", regex="fit_15"),
  check_diagnostics = TRUE,
  format = getOption("cmdstanr_draws_format", NULL)
)

fit_20 <- as_cmdstan_fit(
  fs::dir_ls("Draws/", regex="fit_20"),
  check_diagnostics = TRUE,
  format = getOption("cmdstanr_draws_format", NULL)
)


fit_25 <- as_cmdstan_fit(
  fs::dir_ls("Draws/", regex="fit_25"),
  check_diagnostics = TRUE,
  format = getOption("cmdstanr_draws_format", NULL)
)

# Explore fit ----------------------------------------------------------

fit_15$summary()
fit_20$summary()
fit_25$summary()


draws_15 <- fit_15$draws()
fit15_draws <- as_draws_df(draws_15)
fit15_draws$temp <- 15

draws_20 <- fit_20$draws()
fit20_draws <- as_draws_df(draws_20)
fit20_draws$temp <- 20

draws_25 <- fit_25$draws()
fit25_draws <- as_draws_df(draws_25)
fit25_draws$temp <- 25

FR_posterior_stan <- bind_rows(fit15_draws, fit20_draws, fit25_draws)
FR_posterior_stan$fit <- "Stan_by_temp" 


includepars <- c("b_log", "c_log","h_log", "K_log", "q", "r_log")
bayesplot::mcmc_trace(draws_25, pars = includepars)
bayesplot::mcmc_areas(draws_25, area_method = "equal height", pars = includepars)
bayesplot::mcmc_pairs(draws_25, pars = includepars)


library(ggplot2)
ggplot(FR_posterior_stan, aes(x=temp, y=(q), group=as.factor(temp))) + 
  geom_violin(trim=FALSE, fill="gray")+
  labs(title="Plot of q  by temp", x="temperature", y = "Estimate")+
  geom_boxplot(width=0.1)+
  theme_classic()





