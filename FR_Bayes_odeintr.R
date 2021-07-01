library(BayesianTools)
library("deSolve")
library("bbmle")
library("emdbook")
library(odeintr)
library(here)

#------------------------------------------------------------------------------
# data 1: Vuciv-Pestic et al.
#------------------------------------------------------------------------------

df = read.csv(here("data/data_vucic-pestic_manual.csv"))
n = nrow(df) 
str(df)
Tt <- rep(1,nrow(df))
P <- rep(1,nrow(df))

FR.gen = ' 
dxdt[0] = -b * pow(x[0],1+q) / (1 + b * h * pow(x[0],1+q)) * P; 
'
compile_sys("FR_gen", FR.gen, pars = c("b","h","q","P"), method = "rk54")

eq.odeint.general = function(N0, b, h, q, P, Tt, timesteplength){
  FR_gen_set_params(b = b, h = h, q = q, P = P)
  FR_gen(N0,Tt,timesteplength)
}

eaten.odeint.general = function(N0 = df$N0, b, h, q, Tt = rep(1,nrow(df)), P = rep(1,nrow(df)), steps=100){
  Neaten.est = vector()
  for(i.eaten in 1:length(N0)){
    Neaten.est[i.eaten] = N0[i.eaten] - eq.odeint.general(N0 = N0[i.eaten],
                                                          b = b,
                                                          h = h,
                                                          q = q,
                                                          P = P[i.eaten],
                                                          Tt= Tt[i.eaten],
                                                          timesteplength=Tt[i.eaten]/steps)[steps+1,2]
  }
  return(Neaten.est)
}



nll.odeint.general = function(parms, Neaten = df$Neaten, N0 = df$N0, Tt = rep(1,nrow(df)), P = rep(1,nrow(df)), steps=100){
  #if(b <= 0 || h <= 0) return(Inf)
  y = eaten.odeint.general(N0=N0, b=parms[1], h=parms[2], q=parms[3], Tt=Tt, P=P, steps=steps)

  # Bayesian tools expect likelihood, not NLL. Likelihood function adjusted here
    ll = sum(dbinom(x = Neaten, 
                      size = N0, 
                      prob = y/N0, 
                      log=T))
  return(ll)
}


refPars <- data.frame(best=c(0.5,0.2,1), lower=c(0.001, 0.0001, -.5), upper=c(1, 0.2, 2), row.names=c("b", "h", "q"))
prior <- createUniformPrior(lower = refPars$lower, 
                            upper = refPars$upper, 
                            best = refPars$best)
bayesianSetup <- createBayesianSetup(nll.odeint.general, prior=prior, names=c("b", "h", "q"))
iter = 1000
settings = list(iterations = iter, message = F)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

## Not run: 
plot(out)
summary(out)
marginalPlot(out)
correlationPlot(out)
gelmanDiagnostics(out) # should be below 1.05 for all parameters to demonstrate convergence 









