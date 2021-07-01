#------------------------------------------------------------------------------
# Fitting functional responses - Direct parameter estimation by simulating differential equations
# 
# R-code, source functions for odeintr
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
# generalized FR
#------------------------------------------------------------------------------

FR.gen = ' 
dxdt[0] = -b * pow(x[0],1+q) / (1 + b * h * pow(x[0],1+q)) * P; 
'
compile_sys("FR_gen", FR.gen, pars = c("b","h","q","P"), method = "rk54")

eq.odeint.general = function(N0, b, h, q, P, Tt, timesteplength){
  FR_gen_set_params(b = b, h = h, q = q, P = P)
  FR_gen(N0,Tt,timesteplength)
}

eaten.odeint.general = function(N0, b, h, q, Tt, P, steps=100){
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

nll.odeint.general = function(Neaten, N0, b, h, q, Tt, P, steps=100){
  if(b <= 0 || h <= 0) return(Inf)
  y = eaten.odeint.general(N0=N0, b=b, h=h, q=q, Tt=Tt, P=P, steps=steps)
  nll = -1*sum(dbinom(x = Neaten, 
                      size = N0, 
                      prob = y/N0, 
                      log=T))
  return(nll)
}

