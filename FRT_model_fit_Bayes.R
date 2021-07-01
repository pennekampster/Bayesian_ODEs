# ----------------------------------------------------------------------------------------------------------------------------------
# Functional response estimation
# ----------------------------------------------------------------------------------------------------------------------------------
#
# The following R-code carries out the model fitting done in: https://www.biorxiv.org/content/10.1101/498030v2
#
# The code extends (and is based on):
#
# Rosenbaum, B. & Rall, B.C. (2018). Fitting functional responses: Direct parameter estimation by
# simulating differential equations. Methods in Ecology and Evolution.
#
# ----------------------------------------------------------------------------------------------------------------------------------
#
# Uriah Daugaard
#
# Department of Evolutionary Biology and Environmental Studies, University of Zurich, Switzerland
#
# May 2019
#
# ----------------------------------------------------------------------------------------------------------------------------------

############## PACKAGES
library(odeintr)
library(dplyr)
library(readr)
library(bbmle)

############## DATA
FRT_Dataset <- read_csv(here("data/FRT_Dataset.csv"))
FRT_Dataset_15 <- filter(FRT_Dataset, Temperature==15)
FRT_Dataset_20 <- filter(FRT_Dataset, Temperature==20)
FRT_Dataset_25 <- filter(FRT_Dataset, Temperature==25)

############## FUNCTIONS
source("FRT_function_Bayes.R")

############## FITTING

# | TEMPERATURE  15°C |

fit.15 = mle2(minuslogl = nll.odeint.general.pred,
              start = list(b.log = -4.106669,
                           h = exp(-3.531813),
                           q = 0.867401,
                           r = exp(-0.994875),
                           K.log = 7.901859,
                           c = exp(-5.600194),
                           sigma = 0.222916),
              data = with(list(N0 = Prey_start_density,
                               Ndead = Prey_start_density - Prey_end_density,
                               P = Predator_start_density,
                               P.end = Predator_end_density,
                               Tt = Incubation_time),
                          data = FRT_Dataset_15))
summary(fit.15)


refPars <- data.frame(best=c(-4.106669, exp(-3.531813), 0.867401, exp(-0.994875), 7.901859, exp(-5.600194), 0.222916), 
                      lower = c(-10, 0, -0.99, 0, 1, 0, 0),
                      upper = c(-1, 0.1, 3, 1, 10, 0.1, 0.3),
                      row.names=c("b.log", "h", "q", "r", "K.log", "c", "sigma"))


prior <- createUniformPrior(lower = refPars$lower, 
                            upper = refPars$upper, 
                            best = refPars$best)

bayesianSetup <- createBayesianSetup(nll.odeint.general.pred, prior=prior, names=c("b.log", "h", "q", "r", "K.log", "c", "sigma"))
iter = 1000
settings = list(iterations = iter, message = F)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

## Not run: 
plot(out)
summary(out)
marginalPlot(out)
correlationPlot(out)
gelmanDiagnostics(out) # should be below 1.05 for all parameters to demonstrate convergence 








# | TEMPERATURE  20°C |

fit.20 = mle2(minuslogl = nll.odeint.general.pred,
              start = list(b.log = -15.229665,
                           h = exp(-3.100683), 
                           q = 3.492671, 
                           r = exp(0.136438),
                           K.log = 7.897008, 
                           c = exp(-4.554020), 
                           sigma = 0.236063),
              data = with(list(N0 = Prey_start_density,
                               Ndead = Prey_start_density - Prey_end_density,
                               P = Predator_start_density,
                               P.end = Predator_end_density,
                               Tt = Incubation_time),
                          data = FRT_Dataset_20))
summary(fit.20)


# | TEMPERATURE  25°C |

fit.25 = mle2(minuslogl = nll.odeint.general.pred,
              start = list(b.log = 1.116512, 
                           h = exp(-6.987074),
                           q=-0.569485,
                           r = exp(0.189600),
                           K.log = 8.476514,
                           c = exp(-4.580168),
                           sigma = 0.422386),
              data = with(list(N0 = Prey_start_density,
                               Ndead = Prey_start_density - Prey_end_density,
                               P = Predator_start_density,
                               P.end = Predator_end_density,
                               Tt = Incubation_time),
                          data = FRT_Dataset_25))
summary(fit.25)