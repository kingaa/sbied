## ----opts,include=FALSE--------------------------------------------------
library(pomp)
library(knitr)
prefix <- "polio"
opts_chunk$set(
  progress=TRUE,
  prompt=FALSE,tidy=FALSE,highlight=TRUE,
  strip.white=TRUE,
  warning=FALSE,
  message=FALSE,
  error=FALSE,
  echo=TRUE,
  cache=TRUE,
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

## ----data----------------------------------------------------------------
polio_data <- read.table("polio_wisconsin.csv")
colnames(polio_data)

## ----package-------------------------------------------------------------
require(pomp)
packageVersion("pomp")

## ----statenames----------------------------------------------------------
polio_statenames <- c("SB1","SB2","SB3","SB4","SB5","SB6","IB","SO","IO")
polio_obsnames <- "cases"
polio_t0 <- 1932+4/12

## ----covariates----------------------------------------------------------
polio_K <- 6
polio_tcovar <- polio_data$time
polio_bspline_basis <- periodic.bspline.basis(polio_tcovar,nbasis=polio_K,degree=3,period=1)
colnames(polio_bspline_basis)<- paste("xi",1:polio_K,sep="")
covartable <- data.frame(
  time=polio_tcovar,
  polio_bspline_basis,
  B=polio_data$births,
  P=predict(smooth.spline(x=1931:1954,y=polio_data$pop[12*(1:24)]),
               x=polio_tcovar)$y
)

## ----rp_names------------------------------------------------------------
polio_rp_names <- c("b1","b2","b3","b4","b5","b6","psi","rho","tau","sigma_dem","sigma_env")

## ----ivp_names-----------------------------------------------------------
polio_ivp_names <- c("SO_0","IO_0")
polio_paramnames <- c(polio_rp_names,polio_ivp_names)

## ----fixed_names---------------------------------------------------------
polio_fp_names <- c("delta","K","SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0")
polio_paramnames <- c(polio_rp_names,polio_ivp_names,polio_fp_names)

## ----polio_check_fixed---------------------------------------------------
covar_index_t0 <- which(abs(covartable$time-polio_t0)<0.01)
polio_initial_births <- as.numeric(covartable$B[covar_index_t0-0:5])
names(polio_initial_births) <- c("SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0") 
polio_fixed_params <- c(delta=1/60,K=polio_K,polio_initial_births)

## ----polio_read_mle------------------------------------------------------
polio_params <- data.matrix(read.table("polio_params.csv",row.names=NULL,header=TRUE))
polio_mle <- polio_params[which.max(polio_params[,"logLik"]),][polio_paramnames]
polio_mle[polio_fp_names]
polio_fixed_params

## ----rprocess------------------------------------------------------------
polio_rprocess <- Csnippet("
  double lambda, beta, var_epsilon, p, q;
 
  beta = exp(dot_product( (int) K, &xi1, &b1));
  lambda = (beta * (IO+IB) / P + psi);
  var_epsilon = pow(sigma_dem,2)/ lambda +  pow(sigma_env,2);
  lambda *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon,var_epsilon);
  p = exp(- (delta+lambda)/12);
  q = (1-p)*lambda/(delta+lambda);
  SB1 = B;
  SB2= SB1*p;
  SB3=SB2*p;
  SB4=SB3*p;
  SB5=SB4*p;
  SB6=SB5*p;
  SO= (SB6+SO)*p;
  IB=(SB1+SB2+SB3+SB4+SB5+SB6)*q;
  IO=SO*q;
")

## ----measure-------------------------------------------------------------
polio_dmeasure <- Csnippet("
  double tol = 1.0e-25;
  double mean_cases = rho*IO;
  double sd_cases = sqrt(pow(tau*IO,2) + mean_cases);
  if(cases > 0.0){
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
  } else{
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
  }
  if (give_log) lik = log(lik);
")

polio_rmeasure <- Csnippet("
  cases = rnorm(rho*IO, sqrt( pow(tau*IO,2) + rho*IO ) );
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")

## ----initializer---------------------------------------------------------
polio_initializer <- Csnippet("
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

## ----trans---------------------------------------------------------------
polio_toEstimationScale <- Csnippet("
 Tpsi = log(psi);
 Trho = logit(rho);
 Ttau = log(tau);
 Tsigma_dem = log(sigma_dem);
 Tsigma_env = log(sigma_env);
 TSO_0 =  logit(SO_0);
 TIO_0 = logit(IO_0);
")

polio_fromEstimationScale <- Csnippet("
 Tpsi = exp(psi);
 Trho = expit(rho);
 Ttau = exp(tau);
 Tsigma_dem = exp(sigma_dem);
 Tsigma_env = exp(sigma_env);
 TSO_0 =  expit(SO_0);
 TIO_0 = expit(IO_0);
")

## ----pomp----------------------------------------------------------------
polio <- pomp(
  data=subset(polio_data, select=c("cases","time"),(time > polio_t0 + 0.01) & (time < 1953+1/12+0.01)),
  times="time",
  t0=polio_t0,
  params=polio_mle,
  rprocess = euler.sim(step.fun = polio_rprocess, delta.t=1/12),
  rmeasure= polio_rmeasure,
  dmeasure = polio_dmeasure,
  covar=covartable,
  tcovar="time",
  obsnames = polio_obsnames,
  statenames = polio_statenames,
  paramnames = polio_paramnames,
  covarnames = c("xi1","B","P"),
  initializer=polio_initializer,
  toEstimationScale=polio_toEstimationScale, 
  fromEstimationScale=polio_fromEstimationScale
) 
plot(polio)

## ----run_level-----------------------------------------------------------
run_level <- 1
polio_Np <-          c(100,5e3,1e4)
polio_Nmif <-        c(10, 200,400)
polio_Nreps_eval <-  c(2,  10,  20)
polio_Nreps_local <- c(10, 20, 40)
polio_Nreps_global <-c(10, 20, 100)
polio_Nsim <-        c(50,100, 500) 

## ----pf1-----------------------------------------------------------------
require(doParallel)
cores <- 20
registerDoParallel(cores)
t1 <- system.time(
  pf1 <- foreach(i=1:cores,.packages='pomp') %dopar% try(
    pfilter(polio,Np=polio_Np[run_level],seed=297221+i)
  )
)
(L1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))

## ----persistence---------------------------------------------------------
t_sim <- system.time(
 sim <- foreach(i=1:polio_Nsim[run_level],.packages='pomp') %dopar% simulate(polio,seed=1099999+i)
)

no_cases_data <- sum(obs(polio)==0)
no_cases_sim <- sum(sapply(sim,obs)==0)/length(sim)
fadeout1_sim <- sum(sapply(sim,function(po)states(po)["IB",]+states(po)["IO",]<1))/length(sim)
fadeout100_sim <- sum(sapply(sim,function(po)states(po)["IB",]+states(po)["IO",]<100))/length(sim)
imports_sim <- coef(polio)["psi"]*mean(sapply(sim,function(po) mean(states(po)["SO",]+states(po)["SB1",]+states(po)["SB2",]+states(po)["SB3",]+states(po)["SB4",]+states(po)["SB5",]+states(po)["SB6",])))/12

## ----plot_simulated------------------------------------------------------
mle_simulation <- simulate(polio,seed=127)
plot(mle_simulation)

## ----mif-----------------------------------------------------------------
polio_rw.sd_rp <- 0.02
polio_rw.sd_ivp <- 0.2
polio_cooling.fraction.50 <- 0.5

t2 <- system.time({
 m2 <- foreach(i=1:polio_Nreps_local[run_level],
   .packages='pomp', .combine=c) %dopar% try(
   mif2(polio,
     seed=143275+i,
     Np=polio_Np[run_level],
     Nmif=polio_Nmif[run_level],
     cooling.type="geometric",
     cooling.fraction.50=polio_cooling.fraction.50,
     transform=TRUE,
     rw.sd=rw.sd(
       b1=polio_rw.sd_rp,
       b2=polio_rw.sd_rp,
       b3=polio_rw.sd_rp,
       b4=polio_rw.sd_rp,
       b5=polio_rw.sd_rp,
       b6=polio_rw.sd_rp,
       psi=polio_rw.sd_rp,
       rho=polio_rw.sd_rp,
       tau=polio_rw.sd_rp,
       sigma_dem=polio_rw.sd_rp,
       sigma_env=polio_rw.sd_rp,
       IO_0=ivp(polio_rw.sd_ivp),
       SO_0=ivp(polio_rw.sd_ivp)
     )
   )
 )

 lik_m2 <- foreach(i=1:polio_Nreps_local[run_level],.packages='pomp',.combine=rbind) %dopar% {
   set.seed(87932+i)
   logmeanexp(replicate(polio_Nreps_eval[run_level], logLik(pfilter(polio,params=coef(m2[[i]]),Np=polio_Np[run_level]))),
     se=TRUE)
 }
})

r2 <- data.frame(logLik=lik_m2[,1],logLik_se=lik_m2[,2],t(sapply(m2,coef)))
if(run_level>1) write.table(r2,file="polio_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r2$logLik,digits=5)

## ----pairs---------------------------------------------------------------
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(r2,logLik>max(logLik)-20))

## ----box-----------------------------------------------------------------
polio_box <- rbind(
  b1=c(-2,8),
  b2=c(-2,8),
  b3=c(-2,8),
  b4=c(-2,8),
  b5=c(-2,8),
  b6=c(-2,8),
  psi=c(0,0.1),
  rho=c(0,0.1),
  tau=c(0,0.1),
  sigma_dem=c(0,0.5),
  sigma_env=c(0,1),
  SO_0=c(0,1),
  IO_0=c(0,0.01)
)

## ----box_eval------------------------------------------------------------
t3 <- system.time({
 m3 <- foreach(i=1:polio_Nreps_global[run_level],.packages='pomp', .combine=c) %dopar%  mif2(
      m2[[1]],
      seed=1587690+i, 
      start=c(apply(polio_box,1,function(x)runif(1,x)),polio_fixed_params)
 )
 lik_m3 <- foreach(i=1:polio_Nreps_global[run_level],.packages='pomp',.combine=rbind) %dopar% {
   set.seed(87932+i)
   logmeanexp(replicate(polio_Nreps_eval[run_level], logLik(pfilter(polio,params=coef(m3[[i]]),Np=polio_Np[run_level]))), se=TRUE)
 }
})

r3 <- data.frame(logLik=lik_m3[,1],logLik_se=lik_m3[,2],t(sapply(m3,coef)))
if(run_level>1) write.table(r3,file="polio_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r3$logLik,digits=5)

## ----pairs_global--------------------------------------------------------
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(r3,logLik>max(logLik)-20))

## ----nbinom--------------------------------------------------------------
nb_lik <- function(theta) -sum(dnbinom(as.vector(obs(polio)),size=exp(theta[1]),prob=exp(theta[2]),log=TRUE))
nb_mle <- optim(c(0,-5),nb_lik)
-nb_mle$value

## ----arma----------------------------------------------------------------
log_y <- log(as.vector(obs(polio))+1)
arma_fit <- arima(log_y,order=c(2,0,2),seasonal=list(order=c(1,0,1),period=12))
arma_fit$loglik-sum(log_y)

## ----param_file----------------------------------------------------------
polio_params <- read.table("polio_params.csv",row.names=NULL,header=TRUE)
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(polio_params,logLik>max(logLik)-20))

## ----global_rho----------------------------------------------------------
plot(logLik~rho,data=subset(r3,logLik>max(r3$logLik)-10),log="x")

## ----check_class_m3------------------------------------------------------
class(m3)

## ----check_class_m3_1----------------------------------------------------
class(m3[[1]])

## ----mif_diagnostics-----------------------------------------------------
plot(m3[r3$logLik>max(r3$logLik)-10])

## ----likelihood_convergence----------------------------------------------
loglik_convergence <- do.call(cbind,conv.rec(m3[r3$logLik>max(r3$logLik)-10],"loglik"))
matplot(loglik_convergence,type="l",lty=1,ylim=max(loglik_convergence,na.rm=T)+c(-10,0))

## ----save, include=FALSE-------------------------------------------------
if(run_level>1) save(list = ls(all = TRUE), file = "Rout.rda")

