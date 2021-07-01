library(BayesianTools)
library("deSolve")
library("bbmle")
library("emdbook")


#------------------------------------------------------------------------------
# data 1: Vuciv-Pestic et al.
#------------------------------------------------------------------------------

df = read.csv(here("data/data_vucic-pestic_manual.csv"))
n = nrow(df) 
str(df)
Tt <- rep(1,nrow(df))
P <- rep(1,nrow(df))

eq.ode.general = function(t, x, parms){
  with(as.list(parms),{
    dN = -b * x[1]^(1+q) / (1 + b * h * x[1]^(1+q)) * P
    return(list(c(dN)))
  })
}

# define functions
eaten.ode.general = function(N0 = df$N0, b, h, q, Tt = rep(1,nrow(df)), P = rep(1,nrow(df)), steps=100){
  Neaten.est = vector()
  for(i.eaten in 1:length(N0)){
    Neaten.est[i.eaten] = N0[i.eaten] - lsoda(y = N0[i.eaten],
                                              times = seq(0,Tt[i.eaten],length=steps),
                                              func = eq.ode.general,
                                              parms = c(b=b, h=h, q=q, P=P[i.eaten])
    )[length(seq(0,Tt[i.eaten],length=steps)),2]
  }
  return(Neaten.est)
}

#eaten.ode.general(b=b.gen, h=h.gen, q=q.gen)

nll.ode.general = function(Neaten, N0, b, h, q, Tt, P, steps=100){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.ode.general(N0=N0, b=b, h=h, q=q, Tt=Tt, P=P, steps=steps)
  nll = -1*sum(dbinom(x = Neaten, 
                      size = N0, 
                      prob = y/N0, 
                      log=T))
  return(nll)
}

#nll.ode.general()



fit.gen = mle2(minuslogl = nll.ode.general,
               start = list(b = 1,
                            h = 1/25,
                            q = 0),
               data = list(Neaten = df$Neaten,
                           N0 = df$N0,
                           P = rep(1,n),
                           Tt = rep(1,n))
)

summary(fit.gen)

b.gen = summary(fit.gen)@coef[1]
h.gen = summary(fit.gen)@coef[2]
q.gen = summary(fit.gen)@coef[3]

N0s = seq(0, max(df$N0), 1)

Ne.gen = eaten.ode.general(N0 = N0s, 
                           P = rep(1,length(N0s)), 
                           Tt = rep(1,length(N0s)),
                           b = b.gen, 
                           h = h.gen, 
                           q = q.gen)

plot(df)
lines(N0s, Ne.gen, lty=1)


nll.ode.general2 = function(parms, Neaten = df$Neaten, N0 = df$N0, Tt = rep(1,nrow(df)), P = rep(1,nrow(df)), steps=100){
  #parms = refPars$best
  # if(pars[1] <= 0 || pars[2] <= 0) return(Inf)
  y = eaten.ode.general(N0=N0, b=parms[1], h=parms[2], q=parms[3], Tt=Tt, P=P, steps=steps)
  ll = sum(dbinom(x = Neaten, 
                  size = N0, 
                  prob = y/N0, 
                  log=T))
  return(ll)
}


refPars <- data.frame(best=c(b.gen,h.gen,q.gen), lower=c(0.001, 0.0001, -.5), upper=c(1, 0.2, 2), row.names=c("b", "h", "q"))
prior <- createUniformPrior(lower = refPars$lower, 
                            upper = refPars$upper, 
                            best = refPars$best)
bayesianSetup <- createBayesianSetup(nll.ode.general2, prior=prior, names=c("b", "h", "q"))
iter = 1000
settings = list(iterations = iter, message = F)

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

## Not run: 
plot(out)
summary(out)
marginalPlot(out)
correlationPlot(out)
gelmanDiagnostics(out) # should be below 1.05 for all parameters to demonstrate convergence 









