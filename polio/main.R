library(pomp)
library(tidyverse)
library(doParallel)
library(doRNG)
stopifnot(packageVersion("pomp")>="3.4.5")
options(dplyr.summarise.inform=FALSE)

library(tidyverse)
data <- read_csv(
  "https://kingaa.github.io/sbied/polio/polio_wisconsin.csv",
  comment="#")
head(data,5)



statenames <- c("SB1","SB2","SB3","SB4","SB5","SB6","IB","SO","IO")
t0 <- 1932+4/12

library(pomp)
K <- 6
covar <- covariate_table(
  t=data$time,
  B=data$births,
  P=predict(smooth.spline(x=1931:1954,
    y=data$pop[seq(12,24*12,by=12)]))$y,
  periodic.bspline.basis(t,nbasis=K,
    degree=3,period=1,names="xi%d"),
  times="t"
)

rp_names <- c("b1","b2","b3","b4","b5","b6",
  "psi","rho","tau","sigma_dem","sigma_env")

ivp_names <- c("SO_0","IO_0")
paramnames <- c(rp_names,ivp_names)

fp_names <- c("delta","K",
  "SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0")
paramnames <- c(rp_names,ivp_names,fp_names)
covar_index_t0 <- which(abs(covar@times-t0)<0.01)
initial_births <- covar@table["B",covar_index_t0-0:5]
names(initial_births) <- c("SB1_0","SB2_0",
  "SB3_0","SB4_0","SB5_0","SB6_0") 
fixed_params <- c(delta=1/60,K=K,initial_births)

params_guess <- c(
  b1=3,b2=0,b3=1.5,b4=6,b5=5,b6=3,
  psi=0.002,rho=0.01,tau=0.001,
  sigma_dem=0.04,sigma_env=0.5,
  SO_0=0.12,IO_0=0.001,
  fixed_params)

rprocess <- Csnippet("
double beta = exp(dot_product( (int) K, &xi1, &b1));
double lambda = (beta * (IO+IB) / P + psi);
double var_epsilon = pow(sigma_dem,2)/ lambda +  
  pow(sigma_env,2);
lambda *= (var_epsilon < 1.0e-6) ? 1 : 
  rgamma(1/var_epsilon,var_epsilon);
double p = exp(-(delta+lambda)/12);
double q = (1-p)*lambda/(delta+lambda);
SB1=B;
SB2=SB1*p;
SB3=SB2*p;
SB4=SB3*p;
SB5=SB4*p;
SB6=SB5*p;
SO=(SB6+SO)*p;
IB=(SB1+SB2+SB3+SB4+SB5+SB6)*q;
IO=SO*q;
")

dmeasure <- Csnippet("
double tol = 1.0e-25;
double mean_cases = rho*IO;
double sd_cases = sqrt(pow(tau*IO,2) + mean_cases);
if(cases > 0.0){
  lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0)
    - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
} else{
  lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
}
if (give_log) lik = log(lik);")
rmeasure <- Csnippet("
cases = rnorm(rho*IO, sqrt( pow(tau*IO,2) + rho*IO ) );
if (cases > 0.0) {
  cases = nearbyint(cases);
} else {
  cases = 0.0;
}")

rinit <- Csnippet("
SB1 = SB1_0;
SB2 = SB2_0;
SB3 = SB3_0;
SB4 = SB4_0;
SB5 = SB5_0;
SB6 = SB6_0;
IB = 0;
IO = IO_0 * P;
SO = SO_0 * P;
")

partrans <- parameter_trans(
  log=c("psi","rho","tau","sigma_dem","sigma_env"),
  logit=c("SO_0","IO_0")
)

data %>%
  filter(
    time > t0 + 0.01,
    time < 1953+1/12+0.01
  ) %>%
  select(cases,time) %>%
  pomp(
    times="time",t0=t0,
    params=params_guess,
    rprocess=euler(step.fun=rprocess,delta.t=1/12),
    rmeasure=rmeasure,
    dmeasure=dmeasure,
    rinit=rinit,
    partrans=partrans,
    covar=covar,
    statenames=statenames,
    paramnames=paramnames
  ) -> polio

simulate(polio)

run_level <- 3
Np <-          switch(run_level,100, 1e3, 5e3)
Nmif <-        switch(run_level, 10, 100, 200)
Nreps_eval <-  switch(run_level,  2,  10,  20)
Nreps_local <- switch(run_level, 10,  20,  40)
Nreps_global <-switch(run_level, 10,  20, 100)
Nsim <-        switch(run_level, 50, 100, 500) 

library(doParallel)
registerDoParallel()
library(doRNG)

stew(file="results/pf1.rda",{
  registerDoRNG(3899882)
  pf1 <- foreach(i=1:20,.packages="pomp",
    .export=c("polio","Np")) %dopar%
    pfilter(polio,Np=Np)
},info=TRUE,dependson=Np)
L1 <- logmeanexp(sapply(pf1,logLik),se=TRUE)

simulate(polio,nsim=Nsim,seed=1643079359,
  format="data.frame",include.data=TRUE) -> sims

sims %>%
  group_by(.id) %>%
  summarize(
    no_cases=sum(cases==0),
    fadeout1=sum(IO+IB<1,na.rm=TRUE),
    fadeout100=sum(IO+IB<100,na.rm=TRUE),
    imports=coef(polio,"psi")*mean(SO+SB1+SB2+SB3+SB4+SB5+SB6,na.rm=TRUE)/12
  ) %>%
  ungroup() %>%
  gather(var,val,-.id) %>%
  group_by(
    type=if_else(.id=="data","data","sim"),
    var
  ) %>%
  summarize(val=mean(val)) %>%
  ungroup() %>%
  spread(var,val) -> summ
summ %>%
  column_to_rownames("type") -> summ



rw_sd <- eval(substitute(rw.sd(
  b1=rwr,b2=rwr,b3=rwr,b4=rwr,b5=rwr,b6=rwr,
  psi=rwr,rho=rwr,tau=rwr,sigma_dem=rwr,
  sigma_env=rwr,
  IO_0=ivp(rwi),SO_0=ivp(rwi)),
  list(rwi=0.2,rwr=0.02)))

exl <- c("polio","Np","Nmif","rw_sd",
  "Nreps_local","Nreps_eval")

stew(file="results/mif.rda",{
  m2 <- foreach(i=1:Nreps_local,
    .packages="pomp",.combine=c,.export=exl) %dopar%
    mif2(polio, Np=Np, Nmif=Nmif, rw.sd=rw_sd,
         cooling.fraction.50=0.5)
  lik_m2 <- foreach(m=m2,.packages="pomp",.combine=rbind,
    .export=exl) %dopar%
    logmeanexp(replicate(Nreps_eval,
      logLik(pfilter(m,Np=Np))),se=TRUE)
},dependson=run_level)

coef(m2) %>% melt() %>% spread(parameter,value) %>%
  select(-.id) %>%
  bind_cols(logLik=lik_m2[,1],logLik_se=lik_m2[,2]) -> r2
r2 %>% arrange(-logLik) %>%
  write_csv("params.csv")
summary(r2$logLik,digits=5)



box <- rbind(
  b1=c(-2,8), b2=c(-2,8),
  b3=c(-2,8), b4=c(-2,8),
  b5=c(-2,8), b6=c(-2,8),
  psi=c(0,0.1), rho=c(0,0.1), tau=c(0,0.1),
  sigma_dem=c(0,0.5), sigma_env=c(0,1),
  SO_0=c(0,1), IO_0=c(0,0.01)
)

if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}

bake(file="results/box_eval1.rds",{
  registerDoRNG(833102018)
  foreach(i=1:Nreps_global,.packages="pomp",
    .combine=c) %dopar%
    mif2(m2[[1]],params=c(fixed_params,
      apply(box,1,function(x)runif(1,x[1],x[2]))))
},dependson=run_level) -> m3
bake(file="results/box_eval2.rds",{
  registerDoRNG(71449038)
  foreach(m=m3,.packages="pomp",
    .combine=rbind) %dopar%
    logmeanexp(replicate(Nreps_eval,
      logLik(pfilter(m,Np=Np))),se=TRUE)
},dependson=run_level) -> lik_m3

coef(m3) %>% melt() %>% spread(parameter,value) %>%
  select(-.id) %>%
  bind_cols(logLik=lik_m3[,1],logLik_se=lik_m3[,2]) -> r3
read_csv("params.csv") %>%
  bind_rows(r3) %>%
  arrange(-logLik) %>%
  write_csv("params.csv")
summary(r3$logLik,digits=5)



nb_lik <- function (theta) {
  -sum(dnbinom(as.numeric(obs(polio)),
               size=exp(theta[1]),prob=exp(theta[2]),log=TRUE))
}
nb_mle <- optim(c(0,-5),nb_lik)
-nb_mle$value

log_y <- log(as.vector(obs(polio))+1)
arma_fit <- arima(log_y,order=c(2,0,2),
  seasonal=list(order=c(1,0,1),period=12))
arma_fit$loglik-sum(log_y)

params <- read_csv("params.csv")
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,
  data=subset(params,logLik>max(logLik)-20))





library(tidyverse)
params %>% 
  filter(logLik>max(logLik)-20) %>%
  select(-logLik,-logLik_se) %>%
  gather(variable,value) %>%
  group_by(variable) %>%
  summarize(min=min(value),max=max(value)) %>%
  ungroup() %>%
  column_to_rownames(var="variable") %>%
  t() -> box

profile_pts <-  switch(run_level,  3,  5,  30)
profile_Nreps <- switch(run_level, 2,  3,  10)

idx <- which(colnames(box)!="rho")
profile_design(
  rho=seq(0.01,0.025,length=profile_pts),
  lower=box["min",idx],upper=box["max",idx],
  nprof=profile_Nreps
) -> starts

profile_rw_sd <- eval(substitute(rw.sd(
  rho=0,b1=rwr,b2=rwr,b3=rwr,b4=rwr,b5=rwr,b6=rwr,
  psi=rwr,tau=rwr,sigma_dem=rwr,sigma_env=rwr,
  IO_0=ivp(rwi),SO_0=ivp(rwi)),
  list(rwi=0.2,rwr=0.02)))

bake(file="results/profile_rho.rds",{  
  registerDoRNG(1888257101)
  foreach(start=iter(starts,"row"),.combine=rbind,
    .packages=c("pomp","dplyr")) %dopar% {
    polio %>% mif2(params=start,
      Np=Np,Nmif=ceiling(Nmif/2),
      cooling.fraction.50=0.5,
      rw.sd=profile_rw_sd
    ) %>%
      mif2(Np=Np,Nmif=ceiling(Nmif/2),
        cooling.fraction.50=0.1
      ) -> mf
    replicate(Nreps_eval,
      mf %>% pfilter(Np=Np) %>% logLik()
    ) %>% logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(logLik=ll[1],logLik.se=ll[2])
  }
},dependson=run_level) -> m4



read_csv("params.csv") %>%
  bind_rows(m4) %>%
  arrange(-logLik) %>%
  write_csv("params.csv",append=TRUE)



plot(m3[r3$logLik>max(r3$logLik)-10])

loglik_convergence <- do.call(cbind,
  traces(m3[r3$logLik>max(r3$logLik)-10],"loglik"))
matplot(loglik_convergence,type="l",lty=1,
  ylim=max(loglik_convergence,na.rm=T)+c(-10,0))
