## ----prelims,include=FALSE,purl=TRUE,cache=FALSE-------------------------
options(
  keep.source=TRUE,
  encoding="UTF-8"
)

## ----load_bbs------------------------------------------------------------
bsflu_data <- read.table("http://kingaa.github.io/sbied/data/bsflu_data.txt")

## ----bsflu_names---------------------------------------------------------
bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("Beta","mu_I","rho","mu_R1","mu_R2")

## ----bsflu_obsnames------------------------------------------------------
colnames(bsflu_data)[1:2]

## ----csnippets_bsflu-----------------------------------------------------
bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-6,give_log);
"

bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-6);
  C = rpois(rho*R2);
"

bsflu_rprocess <- "
  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt*mu_I));
  double t3 = rbinom(R1,1-exp(-dt*mu_R1));
  double t4 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= t1;
  I += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
"

bsflu_fromEstimationScale <- "
 TBeta = exp(Beta);
 Tmu_I = exp(mu_I);
 Trho = expit(rho);
"

bsflu_toEstimationScale <- "
 TBeta = log(Beta);
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
require(pomp)
stopifnot(packageVersion("pomp")>="1.2.2")
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
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  initializer=Csnippet(bsflu_initializer)
)
plot(bsflu)

## ----start_params--------------------------------------------------------
params <- c(Beta=0.005,mu_I=2,rho=0.9,mu_R1=1/3,mu_R2=1/2)

## ----init_sim------------------------------------------------------------
y <- simulate(bsflu,params=params,nsim=10,as.data.frame=TRUE)

require(ggplot2)
require(reshape2)

ggplot(data=melt(subset(y,select=c(time,B,C,sim)),
                 id=c("sim","time")),
       mapping=aes(x=time,y=value,group=sim))+
  geom_line()+
  facet_grid(variable~.)+
  theme_classic()

## ----init_pfilter--------------------------------------------------------
pf <- pfilter(bsflu,params=params,Np=1000)
plot(pf)

## ----run_level-----------------------------------------------------------
run_level <- 3
switch(run_level,
       {bsflu_Np=100; bsflu_Nmif=10; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=20000; bsflu_Nmif=100; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=60000; bsflu_Nmif=300; bsflu_Neval=10; bsflu_Nglobal=100; bsflu_Nlocal=20}
)

## ----fixed_params--------------------------------------------------------
fixed_params <- with(bsflu_data,c(mu_R1=1/(sum(B)/512),mu_R2=1/(sum(C)/512)))
est_params <- params[c("Beta","mu_I","rho")]

## ----parallel-setup,cache=FALSE------------------------------------------
require(foreach)
require(doParallel)

cores <- 20
registerDoParallel(cores)

set.seed(2036049659,kind="L'Ecuyer")
mcopts <- list(set.seed=TRUE)

## ----pf------------------------------------------------------------------
stew(file="pf-1.rda",{
  
  t_pf <- system.time(
    pf <- foreach(i=1:20,.packages='pomp',
                  .options.multicore=list(set.seed=TRUE)) %dopar% 
      try(
        pfilter(bsflu,params=c(est_params,fixed_params),Np=1000)
      )
  )
  
},seed=625904618,kind="L'Ecuyer")

(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))

## ----init_csv------------------------------------------------------------
results <- as.data.frame(as.list(c(est_params,fixed_params,loglik=L_pf[1],loglik.se=L_pf[2])))
write.csv(results,file="mif_bsflu_params.csv",row.names=FALSE)

## ----box_search_local----------------------------------------------------
stew(file="local_search-1.rda",{
  
  t_local <- system.time({
    mifs_local <- foreach(i=1:20,
                          .packages='pomp',
                          .combine=c, 
                          .options.multicore=list(set.seed=TRUE)) %dopar%  
      try(
        mif2(
          bsflu,
          start=c(fixed_params,est_params),
          Np=20000,
          Nmif=100,
          cooling.type="geometric",
          cooling.fraction.50=0.5,
          transform=TRUE,
          rw.sd=rw.sd(Beta=0.02,mu_I=0.02,rho=0.02)
        )
      )
  })  
},seed=482947940,kind="L'Ecuyer")


## ----lik_local_eval------------------------------------------------------
stew(file="lik_local-1.rda",{
  
  t_local_eval <- system.time({
    liks_local <- foreach(mf=mifs_local,
                          .packages='pomp',
                          .combine=rbind) %dopar% 
      try({
        evals <- replicate(10, logLik(pfilter(mf,Np=20000)))
        ll <- logmeanexp(evals,se=TRUE)
        c(coef(mf),loglik=ll)
      })
  })
},seed=900242057,kind="L'Ecuyer")

results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=5)

## ----pairs_local---------------------------------------------------------
pairs(~logLik+Beta+mu_I+rho,data=subset(results_local,logLik>max(logLik)-50))

## ----box-----------------------------------------------------------------
bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_I=c(0.5,2),
  rho = c(0.5,1)
)

## ----box_eval------------------------------------------------------------
stew(file=sprintf("box_eval-%d.rda",run_level),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp', .combine=c, .options.multicore=mcopts) %dopar%  mif2(
      mifs_local[[1]],
      start=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),fixed_params)
    )
  })
},seed=1270401374,kind="L'Ecuyer")

## ----lik_global_eval-----------------------------------------------------
stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp',.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu,params=coef(mifs_global[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))
summary(results_global$logLik,digits=5)

## ----save_params,eval=FALSE----------------------------------------------
## if (run_level>1)
##   write.table(rbind(results_local,results_global),
##               file="mif_bsflu_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)

## ----pairs_global--------------------------------------------------------
pairs(~logLik+Beta+mu_I+rho,data=subset(results_global,logLik>max(logLik)-250))

