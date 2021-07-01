#------------------------------------------------------------------------------
# Fitting functional responses - Direct parameter estimation by simulating differential equations
# 
# R-code
# 
# Benjamin Rosenbaum & Björn Rall
# 
# German Centre for Integrative Biodiversity Research (iDiv) Halle-Jena-Leipzig, Deutscher Platz 5, 04103 Leipzig, Germany
# 
# Institute of Biodiversity, Friedrich Schiller University Jena, Dornburger Str. 159, 07743 Jena, Germany
# 
# Corresponding authors: benjamin.rosenbaum@idiv.de & bjoern.rall@idiv.de
#
# September 2017
#------------------------------------------------------------------------------

rm(list=ls())
	
library("bbmle")
library("emdbook")
library("deSolve")

source("source_fr_manual.R")

#------------------------------------------------------------------------------
# data 1: Vuciv-Pestic et al.
#------------------------------------------------------------------------------

df = read.csv("data_vucic-pestic_manual.csv")
n = nrow(df) 
str(df)

#------------------------------------------------------------------------------
# type II fit
#------------------------------------------------------------------------------

fit.II = mle2(minuslogl = nll.bolker,
              start = list(b = 1,
                           h = 1/25),
              data = list(Neaten = df$Neaten,
                          N0 = df$N0,
                          P = rep(1,n),
                          Tt = rep(1,n))
)

summary(fit.II)

b.II = summary(fit.II)@coef[1]
h.II = summary(fit.II)@coef[2]

N0s = seq(0,max(df$N0),1)

Ne.II = eaten.bolker(N0 = N0s,
                     P = rep(1,length(N0s)),
                     Tt = rep(1,length(N0s)),
                     b = b.II, 
                     h = h.II)

plot(df)
lines(N0s, Ne.II, lty=3)


#------------------------------------------------------------------------------
# type III fit
#------------------------------------------------------------------------------

fit.III = mle2(minuslogl = nll.trexler,
               start = list(b = 1,
                            h = 1/25),
               data = list(Neaten = df$Neaten,
                           N0 = df$N0,
                           P = rep(1,n),
                           Tt = rep(1,n))
)

summary(fit.III)

b.III = summary(fit.III)@coef[1]
h.III = summary(fit.III)@coef[2]

N0s = seq(0,max(df$N0),1)

Ne.III = eaten.trexler(N0 = N0s,
                       P = rep(1,length(N0s)),
                       Tt = rep(1,length(N0s)),
                       b = b.III,
                       h = h.III)

plot(df)
lines(N0s, Ne.III, lty=2)

#------------------------------------------------------------------------------
# generalized fit
#------------------------------------------------------------------------------

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

N0s = seq(0,max(df$N0),1)

Ne.gen = eaten.ode.general(N0 = N0s, 
                           P = rep(1,length(N0s)), 
                           Tt = rep(1,length(N0s)),
                           b = b.gen, 
                           h = h.gen, 
                           q = q.gen)

plot(df)
lines(N0s, Ne.gen, lty=1)

#------------------------------------------------------------------------------
# model comparison
#------------------------------------------------------------------------------

AIC(fit.II)
AIC(fit.III)
AIC(fit.gen)

plot(df)
lines(N0s, Ne.II,  lty=3)
lines(N0s, Ne.III, lty=2)
lines(N0s, Ne.gen, lty=1)
legend("topleft",legend=c("type II","type III","gen."),lty=c(3,2,1))

# abline(1/h.II,  0, lty=3)
# abline(1/h.III, 0, lty=2)
# abline(1/h.gen, 0, lty=1)

#------------------------------------------------------------------------------
# data 2: Archer et al., prey mortality
#------------------------------------------------------------------------------

library("bbmle")
library("deSolve")

source("source_fr_manual.R")

df = read.csv("data_archer_manual.csv")
n = nrow(df)
str(df)

#------------------------------------------------------------------------------
# generalized fit, mortality model
#------------------------------------------------------------------------------

fit.mort = mle2(minuslogl = nll.ode.general.mort,
                start = list(b = 1,
                             h = 1/30,
                             q = 0,
                             m = 0.01),
                data = list(Ndead = df$Ndead,
                            N0 = df$N0,
                            P = df$PredNo,
                            Tt = rep(1,n))
)

summary(fit.mort)

b.mort = fit.mort@coef[[1]]
h.mort = fit.mort@coef[[2]]
q.mort = fit.mort@coef[[3]]
m.mort = fit.mort@coef[[4]]

N0s = seq(0,max(df$N0),1)

Ne.mort = eaten.ode.general.mort(N0 = N0s,
                                 P = rep(1,length(N0s)),
                                 Tt = rep(1,length(N0s)),
                                 b = b.mort,
                                 h = h.mort,
                                 q = q.mort,
                                 m = m.mort)

plot(df$N0[df$PredNo>0],df$Ndead[df$PredNo>0],
     xlab="N0",
     ylab="Ndead")
lines(N0s, Ne.mort)

Ne.mort.control = eaten.ode.general.mort(N0 = N0s,
                                         P = rep(0,length(N0s)),
                                         Tt = rep(1,length(N0s)),
                                         b = b.mort,
                                         h = h.mort,
                                         q = q.mort,
                                         m = m.mort)

plot(df$N0[df$PredNo==0],df$Ndead[df$PredNo==0],
     xlab="N0",
     ylab="Ndead")
lines(N0s, Ne.mort.control)

#------------------------------------------------------------------------------
# data 3: Uszko et al., densities
#------------------------------------------------------------------------------
library("bbmle")
library("deSolve")

source("source_fr_manual.R")

df = read.csv("data_uszko_manual_2.csv")
n=nrow(df)
str(df)

#------------------------------------------------------------------------------
# generalized fit, densities
#------------------------------------------------------------------------------

fit.dens = mle2(minuslogl = nll.ode.general.dens,
                start = list(b = 1,
                             h = 10,
                             q = 0,
                             sigma = 1),
                data = list(Neaten = df$Neaten,
                            N0 = df$N0,
                            Tt = df$T,
                            P = df$P)
)

summary(fit.dens)

b.dens = summary(fit.dens)@coef[1]
h.dens = summary(fit.dens)@coef[2]
q.dens = summary(fit.dens)@coef[3]

N0s = seq(0,max(df$N0),length=100)

Ne.dens = eaten.ode.general(N0 = N0s, 
                            P = rep(df$P[1],length(N0s)), 
                            Tt = rep(df$T[1],length(N0s)),
                            b = b.dens, 
                            h = h.dens, 
                            q = q.dens)

plot(df$N0, df$Neaten)
lines(N0s, Ne.dens, lty=1)

#------------------------------------------------------------------------------
# data 4: Fussmann et al., densities and growth
#------------------------------------------------------------------------------

library("bbmle")
library("deSolve")

source("source_fr_manual.R")

df = read.csv("data_fussman_manual.csv")
n=nrow(df)
str(df)

#------------------------------------------------------------------------------
# generalized fit, densities and growth
#------------------------------------------------------------------------------

# fit.control = mle2(minuslogl = nll.ode.general.growth.dens,
#                    start = list(r.log = log(1e-3),
#                                 K.log = log(5e+5),
#                                 sigma = 1),
#                    data = list(N0 = df$N0[df$P==0],
#                                Ndead = df$Ndead[df$P==0],
#                                P = df$P[df$P==0],
#                                Tt = df$T[df$P==0]),
#                    fixed = list(b.log = -18.87,
#                                 h.log = -1.18,
#                                 q = 0.72)
# )

# fit.treatment = mle2(minuslogl = nll.ode.general.growth.dens,
#                      start = list(b.log = log(1e-6),
#                                   h.log = log(1e-1),
#                                   q = 0.0,
#                                   sigma = 0.19),
#                      data = list(N0 = df$N0[df$P>0],
#                                  Ndead = df$Ndead[df$P>0],
#                                  P = df$P[df$P>0],
#                                  Tt = df$T[df$P>0]),
#                      fixed = list(r.log = -8.43,
#                                   K.log = 12.92)
# )

fit.growth = mle2(minuslogl = nll.ode.general.growth.dens,
                  start = list(b.log = -19.63,
                               h.log = -1.22,
                               q = 0.78,
                               r.log = -8.43,
                               K.log = 12.92,
                               sigma = 0.18),
                    data = list(N0 = df$N0,
                                Ndead = df$Ndead,
                                P = df$P,
                                Tt = df$T)
)

summary(fit.growth)

b.growth = exp(summary(fit.growth)@coef[1])
h.growth = exp(summary(fit.growth)@coef[2])
q.growth = summary(fit.growth)@coef[3]
r.growth = exp(summary(fit.growth)@coef[4])
K.growth = exp(summary(fit.growth)@coef[5])

N0s = seq(0, max(df$N0), length=100)

Ne.growth = eaten.ode.general.growth(N0 = N0s,
                                     P  = rep(200, length(N0s)),
                                     Tt = rep(240, length(N0s)),
                                     b = b.growth,
                                     h = h.growth,
                                     q = q.growth,
                                     r = r.growth,
                                     K = K.growth)
plot(df$N0[df$P>0], df$Ndead[df$P>0],
     xlab="N0",
     ylab="Ndead")
lines(N0s, Ne.growth)

Ne.growth.control = eaten.ode.general.growth(N0 = N0s,
                                             P  = rep(0, length(N0s)),
                                             Tt = rep(240, length(N0s)),
                                             b = b.growth,
                                             h = h.growth,
                                             q = q.growth,
                                             r = r.growth,
                                             K = K.growth)

plot(df$N0[df$P==0], df$Ndead[df$P==0],
     xlab="N0",
     ylab="Ndead")
lines(N0s, Ne.growth.control)
abline(0,0,lty=2)

#------------------------------------------------------------------------------
# odeintr, data 1: Vuciv-Pestic et al.
#------------------------------------------------------------------------------

rm(list=ls())

library("odeintr")
library("bbmle")

source("source_fr_odeintr_manual.R")

df = read.csv("data_vucic-pestic_manual.csv")
n = nrow(df) 

fit.gen = mle2(minuslogl = nll.odeint.general,
               start = list(b = 1,
                            h = 1/25,
                            q = 0),
               data = list(Neaten = df$Neaten,
                           N0 = df$N0,
                           P = rep(1,n),
                           Tt = rep(1,n))
)

summary(fit.gen)





