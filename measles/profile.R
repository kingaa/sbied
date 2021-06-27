## source("codes.R")

if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
} else {
  library(doParallel)
  registerDoParallel()
}

set.seed(594709947L)
library(tidyverse)
library(pomp)
stopifnot(packageVersion("pomp")>="2.1")
theme_set(theme_bw())


















dat %>% 
  pomp(t0=with(dat,2*time[1]-time[2]),
    time="time",
    params=theta,
    rprocess=euler(rproc,delta.t=1/365.25),
    rinit=rinit,
    dmeasure=dmeas,
    rmeasure=rmeas,
    partrans=pt,
    covar=covariate_table(covar,times="time"),
    accumvars=c("C","W"),
    statenames=c("S","E","I","R","C","W"),
    paramnames=c("R0","mu","sigma","gamma","alpha","iota",
      "rho","sigmaSE","psi","cohort","amplitude",
      "S_0","E_0","I_0","R_0")
  ) -> m1


estpars <- setdiff(names(theta),c("sigmaSE","mu","alpha","rho","iota"))

theta["alpha"] <- 1

theta.t <- partrans(m1,theta,"toEst")

theta.t.hi <- theta.t.lo <- theta.t
theta.t.lo[estpars] <- theta.t[estpars]-log(2)
theta.t.hi[estpars] <- theta.t[estpars]+log(2)

profile_design(
  sigmaSE=seq(from=log(0.02),to=log(0.2),length=20),
  lower=theta.t.lo,
  upper=theta.t.hi,
  nprof=40
) -> pd

dim(pd)

pd <- as.data.frame(t(partrans(m1,t(pd),"fromEst")))

pairs(~sigmaSE+R0+mu+sigma+gamma+S_0+E_0,data=pd)

bake("sigmaSE-profile1.rds",{

  library(doRNG)
  registerDoRNG(1598260027L)
  
  foreach (
    p=iter(pd,"row"),
    .combine=bind_rows, .errorhandling="remove", .inorder=FALSE
  ) %dopar% {
    
    tic <- Sys.time()
    
    library(pomp)
    
    m1 %>% 
      mif2(
        params=p,
        Nmif = 50, 
        rw.sd = rw.sd(
          R0=0.02,sigma=0.02,gamma=0.02,psi=0.02,cohort=0.02,amplitude=0.02,
          S_0=ivp(0.02),E_0=ivp(0.02),I_0=ivp(0.02),R_0=ivp(0.02)),
        Np = 1000,
        cooling.type = "geometric",
        cooling.fraction.50 = 0.1
      ) %>%
      mif2() -> mf
    
    ## Runs 10 particle filters to assess Monte Carlo error in likelihood
    pf <- replicate(10, pfilter(mf, Np = 2000))
    ll <- sapply(pf,logLik)
    ll <- logmeanexp(ll, se = TRUE)
    
    toc <- Sys.time()
    etime <- toc-tic
    units(etime) <- "hours"
    
    data.frame(
      as.list(coef(mf)),
      loglik = ll[1],
      loglik.se = ll[2],
      etime = as.numeric(etime)
    )
  }
}) %>%
  filter(is.finite(loglik)) -> sigmaSE_prof



sigmaSE_prof %>%
  mutate(sigmaSE=exp(signif(log(sigmaSE),5))) %>%
  group_by(sigmaSE) %>%
  filter(rank(-loglik)<=20) %>%
  ungroup() -> pd

bake("sigmaSE-profile2.rds",{
  
  library(doRNG)
  registerDoRNG(915963734L)
  
  foreach (p=iter(pd,"row"),
    .combine=rbind, .errorhandling="remove", .inorder=FALSE
  ) %dopar% {
    
    tic <- Sys.time()
    
    library(pomp)
    
    m1 %>% 
      mif2(
        params = p,
        Nmif = 50, 
        rw.sd = rw.sd(
          R0=0.02,sigma=0.02,gamma=0.02,psi=0.02,cohort=0.02,amplitude=0.02,
          S_0=ivp(0.02),E_0=ivp(0.02),I_0=ivp(0.02),R_0=ivp(0.02)),
        Np = 5000,
        cooling.fraction.50 = 0.1
      ) %>%
      mif2() -> mf

    pf <- replicate(10, pfilter(mf, Np = 5000))
    ll <- sapply(pf,logLik)
    ll <- logmeanexp(ll, se = TRUE)
    
    toc <- Sys.time()
    etime <- toc-tic
    units(etime) <- "hours"
    
    data.frame(
      as.list(coef(mf)),
      loglik = ll[1],
      loglik.se = ll[2],
      etime = as.numeric(etime))
  }
}) -> sigmaSE_prof
