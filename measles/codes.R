## ----prelims,include=FALSE,cache=FALSE-----------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
)

stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)
library(tidyverse)
library(pomp)

## ----rproc-------------------------------------------------
rproc <- Csnippet("
  double beta, br, seas, foi, dw, births;
  double rate[6], trans[6];
  
  // cohort effect
  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
    br = cohort*birthrate/dt + (1-cohort)*birthrate;
  else
    br = (1.0-cohort)*birthrate;

  // term-time seasonality
  t = (t-floor(t))*365.25;
  if ((t>=7 && t<=100) ||
      (t>=115 && t<=199) ||
      (t>=252 && t<=300) ||
      (t>=308 && t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
  else
      seas = 1.0-amplitude;

  
  // transmission rate
  beta = R0*(gamma+mu)*seas;

  // expected force of infection
  foi = beta*pow(I+iota,alpha)/pop;
  
  // white noise (extrademographic stochasticity)
  dw = rgammawn(sigmaSE,dt);
  
  rate[0] = foi*dw/dt;  // stochastic force of infection
  rate[1] = mu;         // natural S death
  rate[2] = sigma;      // rate of ending of latent stage
  rate[3] = mu;         // natural E death
  rate[4] = gamma;      // recovery
  rate[5] = mu;         // natural I death

  // Poisson births
  births = rpois(br*dt);
  
  // transitions between classes
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,E,&rate[2],dt,&trans[2]);
  reulermultinom(2,I,&rate[4],dt,&trans[4]);

  S += births   - trans[0] - trans[1];
  E += trans[0] - trans[2] - trans[3];
  I += trans[2] - trans[4] - trans[5];
  R = pop - S - E - I;
  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
  C += trans[4];           // true incidence
")
## ----initializer-------------------------------------------------
rinit <- Csnippet("
  double m = pop/(S_0+E_0+I_0+R_0);
  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);
  W = 0;
  C = 0;
")
## ----dmeasure-------------------------------------------------
dmeas <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 1.0e-18;
  if (cases > 0.0) {
    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
  } else {
    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
  }
")
## ----rmeasure-------------------------------------------------
rmeas <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 1.0e-18;
  cases = rnorm(m,sqrt(v)+tol);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")

  
  
  