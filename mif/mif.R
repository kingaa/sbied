## ----prelims,include=FALSE,purl=TRUE,cache=FALSE-------------------------
library(pomp)
options(cores=20,stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="1.4")

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
library(pomp)
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

library(ggplot2)
library(reshape2)

ggplot(data=melt(subset(y,select=c(time,B,C,sim)),
                 id=c("sim","time")),
       mapping=aes(x=time,y=value,group=sim))+
  geom_line()+
  facet_grid(variable~.)+
  theme_classic()

## ----init_pfilter--------------------------------------------------------
pf <- pfilter(bsflu,params=params,Np=1000)
plot(pf)

## ----fixed_params--------------------------------------------------------
(fixed_params <- with(bsflu_data,c(mu_R1=1/(sum(B)/512),mu_R2=1/(sum(C)/512))))

## ----parallel-setup,cache=FALSE------------------------------------------
library(foreach)
library(doParallel)

registerDoParallel()

set.seed(2036049659,kind="L'Ecuyer")

## ----pf,echo=F,eval=T----------------------------------------------------
stew(file="pf.rda",{
  t_pf <- system.time(
    pf <- foreach(i=1:10,.packages='pomp',
                  .options.multicore=list(set.seed=TRUE)) %dopar% 
      try(
        pfilter(bsflu,params=c(Beta=0.01,mu_I=2,rho=0.9,fixed_params),Np=10000)
      )
  )
},seed=625904618,kind="L'Ecuyer")

(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))

## ----init_csv------------------------------------------------------------
results <- as.data.frame(as.list(c(coef(pf[[1]]),loglik=L_pf[1],loglik=L_pf[2])))
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

## ----box_search_local----------------------------------------------------
stew(file="box_search_local.rda",{
  t_local <- system.time({
    mifs_local <- foreach(i=1:20,
                          .packages='pomp',
                          .combine=c, 
                          .options.multicore=list(set.seed=TRUE)) %dopar%  
      try(
        mif2(
          bsflu,
          start=c(Beta=0.01,mu_I=2,rho=0.9,fixed_params),
          Np=2000,
          Nmif=50,
          cooling.type="geometric",
          cooling.fraction.50=0.5,
          transform=TRUE,
          rw.sd=rw.sd(Beta=0.02,mu_I=0.02,rho=0.02)
        )
      )
  })
},seed=482947940,kind="L'Ecuyer")

ggplot(data=melt(conv.rec(mifs_local)),
       aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()

## ----lik_local-----------------------------------------------------------
stew(file="lik_local.rda",{
  t_local <- system.time({
    results_local <- foreach(mf=mifs_local,
                             .packages='pomp',
                             .combine=rbind,
                             .options.multicore=list(set.seed=TRUE)) %dopar% 
    {
        evals <- replicate(10, logLik(pfilter(mf,Np=20000)))
        ll <- logmeanexp(evals,se=TRUE)
        c(coef(mf),loglik=ll[1],loglik=ll[2])
    }
  })
},seed=900242057,kind="L'Ecuyer")
results_local <- as.data.frame(results_local)

## ----pairs_local---------------------------------------------------------
pairs(~loglik+Beta+mu_I+rho,data=results_local)

## ----local_database------------------------------------------------------
results <- rbind(results,results_local[names(results)])
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

## ----box-----------------------------------------------------------------
bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_I=c(0.5,2),
  rho = c(0.5,1)
)

## ----box_search_global---------------------------------------------------
stew(file="box_search_global.rda",{
  t_global <- system.time({
    results_global <- foreach(i=1:300, 
                              .packages='pomp', 
                              .combine=rbind,
                              .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      guess <- apply(bsflu_box,1,function(x)runif(1,x[1],x[2]))
      mf <- mif2(mifs_local[[1]],start=c(guess,fixed_params),
                 cooling.type='geometric')
      mf <- mif2(mf,Nmif=100)
      ll <- replicate(10,logLik(pfilter(mf,Np=20000)))
      ll <- logmeanexp(ll,se=TRUE)
      c(coef(mf),loglik=ll[1],loglik=ll[2])
    }
  })
},seed=1270401374,kind="L'Ecuyer")
results_global <- as.data.frame(results_global)
results <- rbind(results,results_global[names(results)])
write.csv(results,file="bsflu_params.csv",row.names=FALSE)

## ----pairs_global--------------------------------------------------------
pairs(~loglik+Beta+mu_I+rho,data=results_global)

