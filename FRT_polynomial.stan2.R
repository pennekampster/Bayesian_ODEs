# Library -----------------------------------------------------------------
library(tidyverse)
library(here)
# we recommend running this is a fresh R session or restarting your current session
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
library(readxl)
library(broom)
library(fs)
## make sure to:
#cmdstanr::install_cmdstan(overwrite=T)


# Load and prepare data ---------------------------------------------------------------

D <- read_csv(here("data/FRT_Dataset.csv"))

## returns a list with data structures as declared in "FRT.stan"
getData <- function(D, temp = NULL) {
  if(!is.null(temp)) {
    D <- dplyr::filter(D, Temperature == temp)
  }
  data <- list(
    N_series = nrow(D),
    raw_temp = as.matrix((D["Temperature"]-20)/10),
  #  raw_temp = (as.integer(as.factor(D$Temperature))-2)*0.5,
    time_end = as.matrix(D["Incubation_time"]), # expects array[N_series, 1]
    x_start = cbind(round(D$Prey_start_density), round(D$Predator_start_density)), # expects array[N_series] vector[2]
    x_start_vector = cbind(D$Prey_start_density, D$Predator_start_density), # expects array[N_series] vector[2]
    x_end = cbind(round(D$Prey_end_density), round(D$Predator_end_density))
  )
  return(data)
}


# Fit stan model ----------------------------------------------------------

# get good starting values
# use the estimates of Daugaard et al to fit polynomial and use parameters as starting point for fitting across temperature
parm_temps <- read_excel(here("data", "log_params_temp.xlsx"))
ggplot(data=parm_temps, aes(x=(temperature-20)/10,y=estimate)) + geom_point() + facet_wrap(~term, scales="free_y") + stat_smooth(method="lm", formula = " y ~ poly(x, 2)")
d <-parm_temps %>% mutate(temperature = (temperature-20)/10, temperature_sq = temperature^2) %>% group_by(term) %>% nest() %>% mutate(mod = map(data, ~ lm(estimate ~ temperature  + temperature_sq, data=.)))

# values to initialize the prior
bind_rows(map(d$mod, tidy))


dd <- getData(D, temp = NULL)

model <- cmdstan_model("FRT_polynomial2.stan")

n_chains <- 3
# if(!dir.exists("Draws")) dir.create("Draws")
fit <- model$sample(data = getData(D, temp = NULL),
                       init = replicate(n_chains, list(b2_log = as.numeric(coef(d$mod[[1]])[3]),
                                                       b1_log = as.numeric(coef(d$mod[[1]])[2]),
                                                       b0_log = as.numeric(coef(d$mod[[1]])[1]),
                                                       h2_log = as.numeric(coef(d$mod[[2]])[3]), 
                                                       h1_log = as.numeric(coef(d$mod[[2]])[2]), 
                                                       h0_log = as.numeric(coef(d$mod[[2]])[1]), 
                                                       q2 = as.numeric(coef(d$mod[[3]])[3]),
                                                       q1 = as.numeric(coef(d$mod[[3]])[2]),
                                                       q0 = as.numeric(coef(d$mod[[3]])[1]),
                                                       K2_log = as.numeric(coef(d$mod[[5]])[3]),
                                                       K1_log = as.numeric(coef(d$mod[[5]])[2]),
                                                       K0_log = as.numeric(coef(d$mod[[5]])[1]),
                                                       r2_log = as.numeric(coef(d$mod[[4]])[3]),
                                                       r1_log = as.numeric(coef(d$mod[[4]])[2]),
                                                       r0_log = as.numeric(coef(d$mod[[4]])[1]),
                                                       c2_log = as.numeric(coef(d$mod[[6]])[3]), 
                                                       c1_log = as.numeric(coef(d$mod[[6]])[2]), 
                                                       c0_log = as.numeric(coef(d$mod[[6]])[1]) 
                                                       ), simplify = F), # sigma = 0.222916
                       iter_warmup = 300, iter_sampling = 1000, chains = n_chains, parallel_chains = n_chains, output_dir = "Draws", output_basename = "fit_poly", seed = 1)

fit$save_output_files(dir = "Draws", basename = "fit_poly_one_function_bdf")
#fit <- cmdstanr::read_cmdstan_csv(fs::dir_ls(here("Draws/"), regex="across"))

fit <- as_cmdstan_fit(
  dir_ls("Draws/", regex="fit_poly_one_function"),
  check_diagnostics = TRUE,
  format = getOption("cmdstanr_draws_format", NULL)
)



# fit$save_output_files(dir = "Draws", basename = "fit_25")
# fit <- cmdstanr::read_cmdstan_csv(fit_25$output_files())

# Explore fit ----------------------------------------------------------
# fit <- fit_variational
poly_summary <- fit$summary()

includepars <- c("b_log", "c_log","h_log", "K_log", "q", "r_log")

draws <- fit$draws()

bayesplot::mcmc_trace(draws,regex_pars = "log")
bayesplot::mcmc_areas(draws, area_method = "equal height", regex_pars = "_log")

bayesplot::mcmc_pairs(draws,regex_pars = c("_log","q"),
           off_diag_args = list(size = 1, alpha = 0.5))

fit_draws <- posterior::as_draws_df(draws)


compare  <- cbind(arrange(poly_summary[2:19,c(1,3)], variable), bind_rows(map(d$mod, tidy)))


out <- poly_summary[2:19,c(1,3)]
out$parm <- str_sub(out$variable, 1,1)
out$term <- str_sub(out$variable, 2,2)

pred_df <- pivot_wider(out[,c(2:4)], names_from = "term", values_from = "median")

pred_df <- pred_df %>% rowwise(parm) %>% mutate(pred15 = `2`*(-0.5)^2 + `1` * -0.5 + `0`,
                                     pred20 = `2`*(0)^2 + `1` * 0 + `0`,
                                     pred25 = `2`*(0.5)^2 + `1` * 0.5 + `0`)

preds <- pred_df %>% select(parm, pred15, pred20, pred25) %>% pivot_longer(cols=2:4, names_to="temperature") %>%
  mutate(temperature = as.numeric(gsub("pred", "", temperature)),
         term = ifelse(parm == "q", parm, paste0("log_", parm)))

ggplot() + geom_point(data=parm_temps, aes(x=(temperature-20)/10,y=estimate)) + 
  stat_smooth(data=parm_temps, aes(x=(temperature-20)/10,y=estimate), formula = "y ~ poly(x,2)", method="lm", colur="black")+
  geom_point(data=preds, aes(x=(temperature-20)/10,y=value), colour="red") +
  stat_smooth(data=preds, aes(x=(temperature-20)/10,y=value), colour="red", formula = "y ~ poly(x,2)", method="lm", colur="black")+
  facet_wrap(~term, scales="free_y") 

