#------------------------------------------------------------------------------
# Fitting functional responses - Direct parameter estimation by simulating differential equations
# 
# R-code, source functions
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

#------------------------------------------------------------------------------
# type II FR
#------------------------------------------------------------------------------

eaten.bolker = function(N0, b, h, q=0, Tt, P) {
  a = b*N0^q
  Neaten.est = N0 - lambertW(a*h*N0*exp(-a*(P*Tt-h*N0)))/(a*h)
  return(Neaten.est)
}

nll.bolker = function(Neaten, N0, b, h, q=0, Tt, P){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.bolker(N0=N0, b=b, h=h, q=q, Tt=Tt, P=P)
  nll = -1*sum(dbinom(x = Neaten,
                      size = N0, 
                      prob = y/N0,
                      log=T))
  return(nll)
}

#------------------------------------------------------------------------------
# type III FR
#------------------------------------------------------------------------------

eaten.trexler = function(N0, b, h, Tt, P){
  p = -(b*h*N0^2 + b*P*Tt*N0 + 1) / (b*h*N0)
  q = P*Tt*N0/h
  Neaten.est = -(p/2) - sqrt( (p/2)^2 - q )
  return(Neaten.est)
}

nll.trexler = function(Neaten, N0, b, h, Tt, P){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.trexler(N0=N0, b=b, h=h, Tt=Tt, P=P)
  nll = -1*sum(dbinom(x = Neaten, 
                      size = N0, 
                      prob = y/N0, 
                      log=T))
  return(nll)
}

#------------------------------------------------------------------------------
# generalized FR
#------------------------------------------------------------------------------

eq.ode.general = function(t, x, parms){
  with(as.list(parms),{
    dN = -b * x[1]^(1+q) / (1 + b * h * x[1]^(1+q)) * P
    return(list(c(dN)))
  })
}

eaten.ode.general = function(N0, b, h, q, Tt, P, steps=100){
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

nll.ode.general = function(Neaten, N0, b, h, q, Tt, P, steps=100){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.ode.general(N0=N0, b=b, h=h, q=q, Tt=Tt, P=P, steps=steps)
  nll = -1*sum(dbinom(x = Neaten, 
                      size = N0, 
                      prob = y/N0, 
                      log=T))
  return(nll)
}

#------------------------------------------------------------------------------
# generalized FR and mortality model
#------------------------------------------------------------------------------

eq.ode.general.mort = function(t, x, parms){
  with(as.list(parms),{
    dN = -b * x[1]^(1+q) / (1 + b * h * x[1]^(1+q)) * P - m*x[1]
    return(list(c(dN)))
  })
}

eaten.ode.general.mort = function(N0, b, h, q, m, Tt, P, steps=100){
  Ndead.est = vector()
  for(i.eaten in 1:length(N0)){
    Ndead.est[i.eaten] = N0[i.eaten] - lsoda(y = N0[i.eaten],
                                              times = seq(0,Tt[i.eaten],length=steps),
                                              func = eq.ode.general.mort,
                                              parms = c(b=b, h=h, q=q, m=m, P=P[i.eaten])
    )[length(seq(0,Tt[i.eaten],length=steps)),2]
  }
  return(Ndead.est)
}

nll.ode.general.mort = function(Ndead, N0, b, h, q, m, Tt, P, steps=100){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.ode.general.mort(N0=N0, b=b, h=h, q=q, m=m, Tt=Tt, P=P, steps=steps)
  nll = -1*sum(dbinom(x = Ndead,
                      size = N0,
                      prob = y/N0,
                      log=T))
  return(nll)
}

#------------------------------------------------------------------------------
# generalized FR for densities
#------------------------------------------------------------------------------

nll.ode.general.dens = function(Neaten, N0, b, h, q, Tt, P, steps=100, sigma){
  if(b <= 0 || h <= 0 || sigma <= 0) return(Inf)
  y = eaten.ode.general(N0=N0, b=b, h=h, q=q, Tt=Tt, P=P, steps=steps)
  nll = -1*sum(dnorm(x = log(Neaten), 
                     mean = log(y), 
                     sd = sigma, 
                     log=T))
  return(nll)
}

#------------------------------------------------------------------------------
# generalized FR and growth model for densities
#------------------------------------------------------------------------------

eq.ode.general.growth = function(t, x, parms){
  with(as.list(parms),{
    dN = -b*x[1]^(1+q) / (1+b*h*x[1]^(1+q)) * P + r*x[1]*(1.0-x[1]/K)
    return(list(c(dN)))
  })
}

eaten.ode.general.growth = function(N0, b, h, q, Tt, P, r, K, steps=100){
  eaten.est = vector()
  for(i.eaten in 1:length(N0)){
    eaten.est[i.eaten] = N0[i.eaten] - lsoda(y = N0[i.eaten],
                                             times = seq(0,Tt[i.eaten],length=steps),
                                             func = eq.ode.general.growth,
                                             parms = c(b = b, h = h, q = q, r = r, K = K, P = P[i.eaten])
    )[length(seq(0,Tt[i.eaten],length=steps)),2]
  }
  return(eaten.est)
}

nll.ode.general.growth.dens = function(Ndead, N0, b.log, h.log, q, r.log, K.log, Tt, P, steps=100, sigma){
  if(sigma <= 0 || q <= -1) return(Inf)
  y = eaten.ode.general.growth(N0=N0, 
                               b=exp(b.log), 
                               h=exp(h.log), 
                               q=q, 
                               r=exp(r.log), 
                               K=exp(K.log), 
                               Tt=Tt, 
                               P=P, 
                               steps=steps)
  nll = -1*sum(dnorm(x = log(N0-Ndead), 
                     mean = log(N0-y), 
                     sd = sigma, 
                     log=T))
  return(nll)
}

