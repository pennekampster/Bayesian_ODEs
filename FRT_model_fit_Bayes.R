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
library(here)

############## DATA
FRT_Dataset <- read_csv(here("data/FRT_Dataset.csv"))
FRT_Dataset_15 <- filter(FRT_Dataset, Temperature==15)
FRT_Dataset_20 <- filter(FRT_Dataset, Temperature==20)
FRT_Dataset_25 <- filter(FRT_Dataset, Temperature==25)



############## FUNCTIONS
source("FRT_function_Bayes.R")

############## FITTING

# | TEMPERATURE  15°C |
FRT_Dataset <- FRT_Dataset_15
refPars <- data.frame(best=c(-4.106669, -3.531813,  0.867401, -0.994875, 7.901859, -5.600194,  0.222916), 
                      lower = c(-20, -10, -0.99,  -10, 1,    -10, 0),
                      upper = c(5, 1,     3,   1, 10,    10, 1),
                      row.names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))


# refPars <- data.frame(best=c(-4.106669, exp(-3.531813), 0.867401, exp(-0.994875), 7.901859, -5.600194, 0.222916), 
#                       lower = c(-10, 0, -0.99, 0, 1, 0, 0),
#                       upper = c(-1, 0.1, 3, 1, 10, 0.1, 0.3),
#                       row.names=c("b.log", "h", "q", "r", "K.log", "c", "sigma"))


prior <- createUniformPrior(lower = refPars$lower, 
                            upper = refPars$upper, 
                            best = refPars$best)

#prior <- createTruncatedNormalPrior(refPars$best, rep(2,7), refPars$lower, refPars$upper)

bayesianSetup <- createBayesianSetup(nll.odeint.general.pred, prior=prior, names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))
iter = 1000
settings = list(iterations = iter, message = F)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

# create prior from successful previous run
newPrior = createPriorDensity(out, method = "multivariate",
                              eps = 1e-10, lower = refPars$lower,
                              upper = refPars$upper, best = refPars$best)

bayesianSetup <- createBayesianSetup(nll.odeint.general.pred, prior=newPrior, names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))
iter = 10000
settings = list(iterations = iter, message = F)

out_15 <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

## Not run: 
plot(out_15)
summary(out_15)
marginalPlot(out_15)
correlationPlot(out_15)
gelmanDiagnostics(out_15) # should be below 1.05 for all parameters to demonstrate convergence 


# | TEMPERATURE  20°C |

FRT_Dataset <- FRT_Dataset_20

refPars <- data.frame(best=c(-15.229665, -3.100683,  1.250667, 0.136438, 7.897008, -4.554020,  0.236063), 
                      lower = c(-20, -10, -0.99,  -10, 1,    -10, 0),
                      upper = c(5, 1,     3,   1, 10,    10, 1),
                      row.names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))

# refPars <- data.frame(best=c(-15.229665, exp(-3.100683),  log(3.492671), exp(0.136438), 7.897008, exp(-4.554020),  0.236063), 
#                       lower = c(-10, 0, -0.99, 0, 1, 0, 0),
#                       upper = c(-1, 0.1, 3, 1, 10, 0.1, 0.3),
#                       row.names=c("b.log", "h", "q", "r", "K.log", "c", "sigma"))


prior <- createUniformPrior(lower = refPars$lower, 
                            upper = refPars$upper, 
                            best = refPars$best)

bayesianSetup <- createBayesianSetup(nll.odeint.general.pred, prior=prior, names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))
iter = 1000
settings = list(iterations = iter, message = F)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

# create prior from successful previous run
newPrior = createPriorDensity(out, method = "multivariate",
                              eps = 1e-10, lower = refPars$lower,
                              upper = refPars$upper, best = refPars$best)

bayesianSetup <- createBayesianSetup(nll.odeint.general.pred, prior=newPrior, names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))
iter = 10000
settings = list(iterations = iter, message = F)

out_20 <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)


## Not run: 
plot(out_20)
summary(out_20)
marginalPlot(out_20)
correlationPlot(out_20)
gelmanDiagnostics(out_20) # should be below 1.05 for all parameters to demonstrate convergence 


# | TEMPERATURE  25°C |

# b_log = -4.106669, c_log = -5.600194, h_log = -3.531813, K_log = 7.901859, q = 0.867401, r_log = -0.994875


# refPars <- data.frame(best=c(1, exp(-6.987074),  -0.569485, exp(0.189600), 8.476514, exp(-4.580168),  0.422386), 
#                       lower = c(-20, 0, -0.99, 0, 1,    0, 0),
#                       upper = c(5, 0.1,   3,   2, 10, 0.1, 1),
#                       row.names=c("b.log", "h", "q", "r", "K.log", "c", "sigma"))


FRT_Dataset <- FRT_Dataset_25

refPars <- data.frame(best=c(1, -6.987074,  0.86740, -1.662839, 8.476514, -4.580168,  0.42), 
                      lower = c(-20, -10, -0.99,  -10, 1,    -10, 0),
                      upper = c(5, 1,     3,   1,     10,    10, 1),
                      row.names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))


prior <- createUniformPrior(lower = refPars$lower, 
                            upper = refPars$upper, 
                            best = refPars$best)

bayesianSetup <- createBayesianSetup(nll.odeint.general.pred, prior=prior, names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))
iter = 1000
settings = list(iterations = iter, message = F)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

# create prior from successful previous run
newPrior = createPriorDensity(out, method = "multivariate",
                              eps = 1e-10, lower = refPars$lower,
                              upper = refPars$upper, best = refPars$best)

bayesianSetup <- createBayesianSetup(nll.odeint.general.pred, prior=newPrior, names=c("b.log", "h.log", "q", "r.log", "K.log", "c.log", "sigma"))
iter = 10000
settings = list(iterations = iter, message = F)

out_25 <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)



## Not run: 
plot(out_25)
summary(out_25)
marginalPlot(out_25)
correlationPlot(out_25)
gelmanDiagnostics(out_25) # should be below 1.05 for all parameters to demonstrate convergence 

