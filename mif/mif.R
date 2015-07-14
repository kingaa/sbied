## ----opts,include=FALSE--------------------------------------------------
library(pomp)
library(knitr)
prefix <- "mif"
opts_chunk$set(
  progress=TRUE,
  prompt=FALSE,tidy=FALSE,highlight=TRUE,
  strip.white=TRUE,
  warning=FALSE,
  message=FALSE,
  error=FALSE,
  echo=TRUE,
#  cache=TRUE,
  results='markup',
  fig.show='asis',
  size='small',
  fig.lp="fig:",
  fig.path=paste0("figure/",prefix,"-"),
  cache.path=paste0("cache/",prefix,"-"),
  fig.pos="h!",
  fig.align='center',
  fig.height=4,fig.width=6.83,
  dpi=300,
  dev='png',
  dev.args=list(bg='transparent')
  )

options(
  pomp.cache="cache",
  keep.source=TRUE,
  encoding="UTF-8"
  )

library(ggplot2)
theme_set(theme_bw())

## ----load_bbs------------------------------------------------------------
bsflu_data <- read.table("http://kinglab.eeb.lsa.umich.edu/SBIED/data/bsflu_data.txt")

## ----bsflu_names---------------------------------------------------------
bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("beta","mu_I","rho","mu_R1","mu_R2")

## ----bsflu_obsnames------------------------------------------------------
(bsflu_obsnames <- colnames(bsflu_data)[1:2])

## ----csnippets_bsflu-----------------------------------------------------
bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-6,give_log);
"

bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-6);
  C = rpois(rho*R2);
"

bsflu_rprocess <- "
  double t1 = rbinom(S,1-exp(-beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt*mu_I));
  double t3 = rbinom(R1,1-exp(-dt*mu_R1));
  double t4 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= t1;
  I += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
"

bsflu_fromEstimationScale <- "
 Tbeta = exp(beta);
 Tmu_I = exp(mu_I);
 Trho = expit(rho);
"

bsflu_toEstimationScale <- "
 Tbeta = log(beta);
 Tmu_I = log(mu_I);
 Trho = logit(rho);
"

bsflu_initializer <- "
 S=762;
 I=1;
 R1=0;
 R2=0;
"

## ----pomp_bsflu----------------------------------------------------------
bsflu <- pomp(
     data=bsflu_data,
     times="day",
     t0=0,
     rprocess=euler.sim(
       step.fun=Csnippet(bsflu_rprocess),
       delta.t=1/12
       ),
     rmeasure=Csnippet(bsflu_rmeasure),
     dmeasure=Csnippet(bsflu_dmeasure),
     fromEstimationScale=Csnippet(bsflu_fromEstimationScale),
     toEstimationScale=Csnippet(bsflu_toEstimationScale),
     obsnames = bsflu_obsnames,
     statenames=bsflu_statenames,
     paramnames=bsflu_paramnames,
     initializer=Csnippet(bsflu_initializer)
) 
plot(bsflu)

## ----run_level-----------------------------------------------------------
runlevel <- 3
switch(runlevel,
       {bsflu_Np=100; bsflu_Nmif=10; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=20000; bsflu_Nmif=100; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=60000; bsflu_Nmif=300; bsflu_Neval=10; bsflu_Nglobal=100; bsflu_Nlocal=20}
)

## ----bsflu_params--------------------------------------------------------
bsflu_params <- data.matrix(read.table("mif_bsflu_params.csv",row.names=NULL,header=TRUE))
bsflu_mle <- bsflu_params[which.max(bsflu_params[,"logLik"]),][bsflu_paramnames]

## ----fixed_params--------------------------------------------------------
bsflu_fixed_params <- c(mu_R1=1/(sum(bsflu_data$B)/512),mu_R2=1/(sum(bsflu_data$C)/512))

## ----pf1-----------------------------------------------------------------
require(doParallel)
cores <- 20
registerDoParallel(cores)

## ----pf------------------------------------------------------------------
t_pf <- system.time(
  pf <- foreach(i=1:cores,.packages='pomp') %dopar% try(
    pfilter(bsflu,params=bsflu_mle,Np=bsflu_Np,seed=297221+i)
  )
)
(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))

## ----box_search_local----------------------------------------------------
bsflu_rw.sd <- 0.02
bsflu_cooling.fraction.50 <- 0.5

t_local <- system.time({
  mifs_local <- foreach(i=1:bsflu_Nlocal,.packages='pomp', .combine=c) %dopar%  {
     mif2(
        bsflu,
        start=bsflu_mle,
        seed=743275+i,
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        cooling.type="geometric",
        cooling.fraction.50=bsflu_cooling.fraction.50,
        transform=TRUE,
        rw.sd=rw.sd(
           beta=bsflu_rw.sd,
           mu_I=bsflu_rw.sd,
           rho=bsflu_rw.sd
        )
     )

  }
})

## ----lik_local_eval------------------------------------------------------
t_local_eval <- system.time({
  liks_local <- foreach(i=1:bsflu_Nlocal,.packages='pomp',.combine=rbind) %dopar% {
    set.seed(87932+i)
    evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu,params=coef(mifs_local[[i]]),Np=bsflu_Np)))
    logmeanexp(evals, se=TRUE)
  }
})

results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=5)

## ----pairs_local---------------------------------------------------------
pairs(~logLik+beta+mu_I+rho,data=subset(results_local,logLik>max(logLik)-50))

## ----box-----------------------------------------------------------------
bsflu_box <- rbind(
  beta=c(0.001,0.01),
  mu_I=c(0.5,2),
  rho = c(0.5,1)
)

## ----box_eval------------------------------------------------------------
t_global <- system.time({
  mifs_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp', .combine=c) %dopar%  mif2(
     mifs_local[[1]],
     start=c(apply(bsflu_box,1,function(x)runif(1,x)),bsflu_fixed_params),
     seed=743275+i
  )
})

## ----lik_global_eval-----------------------------------------------------
t_global_eval <- system.time({
  liks_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp',.combine=rbind) %dopar% {
    set.seed(87932+i)
    evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu,params=coef(mifs_global[[i]]),Np=bsflu_Np)))
    logmeanexp(evals, se=TRUE)
  }
})

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)

## ----save_params---------------------------------------------------------
if(runlevel>1) write.table(rbind(results_local,results_global),
   file="mif_bsflu_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)

## ----pairs_global--------------------------------------------------------
pairs(~logLik+beta+mu_I+rho,data=subset(results_global,logLik>max(logLik)-250))

## ----save, include=FALSE-------------------------------------------------
if(runlevel>1) save(list = ls(all = TRUE), file = "Rout.rda")

