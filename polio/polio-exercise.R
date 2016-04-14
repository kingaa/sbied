##' Let's carry out a likelihood evaluation at the reported MLE.
##' To develop and debug code, it is nice to have a version that runs extra quickly, for which we set `run_level=1`.
##' Here, `Np` is the number of particles (i.e., sequential Monte Carlo sample size), and `Nmif` is the number of iterations of the optimization procedure carried out below.
##' Empirically, `Np=5000` and `Nmif=200` are around the minimum required to get stable results with an error in the likelihood of order 1 log unit for this example;
##' this is implemented by setting `run_level=2`.
##' One can then ramp up to larger values for more refined computations, implemented here by `run_level=3`.

set.seed(5996485L)

## ----data----------------------------------------------------------------
polio_data <- read.table("http://kingaa.github.io/sbied/polio/polio_wisconsin.csv")
head(polio_data)

## ----package-------------------------------------------------------------
library(pomp)
packageVersion("pomp")

## ----statenames----------------------------------------------------------
statenames <- c("SB1","SB2","SB3","SB4","SB5","SB6","IB","SO","IO")
t0 <- 1932+4/12

## ----covariates----------------------------------------------------------
bspline_basis <- periodic.bspline.basis(
  polio_data$time,nbasis=6,degree=3,period=1,names="xi%d")
covartable <- data.frame(
  time=polio_data$time,
  B=polio_data$births,
  P=predict(smooth.spline(x=1931:1954,y=polio_data$pop[12*(1:24)]),
            x=polio_data$time)$y,
  bspline_basis
)
head(covartable)

## ----rp_names------------------------------------------------------------
rp_names <- c("b1","b2","b3","b4","b5","b6","psi","rho","tau","sigma_dem","sigma_env")

## ----ivp_names-----------------------------------------------------------
ivp_names <- c("SO_0","IO_0")

## ----fixed_params--------------------------------------------------------
i <- which(abs(covartable$time-t0)<0.01)
initial_births <- as.numeric(covartable$B[i-0:5])
names(initial_births) <- c("SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0") 
fixed_params <- c(delta=1/60,initial_births)
fp_names <- c("delta","SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0")

## ----param_guess---------------------------------------------------------
params <- c(b1=3,b2=0,b3=1.5,b4=6,b5=5,b6=3,psi=0.002,rho=0.01,tau=0.001,
            sigma_dem=0.04,sigma_env=0.5,SO_0=0.12,IO_0=0.001,fixed_params)

## ----rprocess------------------------------------------------------------
rproc <- Csnippet("
  double beta = exp(dot_product(K, &xi1, &b1));
  double lambda = (beta * (IO+IB) / P + psi);
  double var_epsilon = pow(sigma_dem,2)/lambda +  sigma_env*sigma_env;
  lambda *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon,var_epsilon);
  double p = exp(-(delta+lambda)/12);
  double q = (1-p)*lambda/(delta+lambda);
  SB1 = B;
  SB2 = SB1*p;
  SB3 = SB2*p;
  SB4 = SB3*p;
  SB5 = SB4*p;
  SB6 = SB5*p;
  SO = (SB6+SO)*p;
  IB = (SB1+SB2+SB3+SB4+SB5+SB6)*q;
  IO = SO*q;
")

## ----measure-------------------------------------------------------------
dmeas <- Csnippet("
  double tol = 1.0e-25;
  double mean_cases = rho*IO;
  double sd_cases = sqrt(pow(tau*IO,2) + mean_cases);
  if (cases > 0.0) {
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
  } else{
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
  }
  if (give_log) lik = log(lik);
")

rmeas <- Csnippet("
  cases = rnorm(rho*IO, sqrt( pow(tau*IO,2) + rho*IO ) );
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")

## ----initializer---------------------------------------------------------
init <- Csnippet("
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
toEst <- Csnippet("
 Tpsi = log(psi);
 Trho = logit(rho);
 Ttau = log(tau);
 Tsigma_dem = log(sigma_dem);
 Tsigma_env = log(sigma_env);
 TSO_0 =  logit(SO_0);
 TIO_0 = logit(IO_0);
")

fromEst <- Csnippet("
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
  data=subset(polio_data, 
              (time > t0 + 0.01) & (time < 1953+1/12+0.01),	
              select=c("cases","time")),
  times="time",
  t0=t0,
  params=params,
  rprocess = euler.sim(step.fun = rproc, delta.t=1/12),
  rmeasure = rmeas,
  dmeasure = dmeas,
  covar=covartable,
  tcovar="time",
  statenames = statenames,
  paramnames = c(rp_names,ivp_names,fp_names),
  initializer=init,
  toEstimationScale=toEst, 
  fromEstimationScale=fromEst,
  globals="int K = 6;"
)

run_level <- 3
polio_Np <-          c(100,5e3,1e4)
polio_Nmif <-        c(10, 200,400)
polio_Nreps_eval <-  c(2,  10,  20)
polio_Nreps_local <- c(10, 20, 40)
polio_Nreps_global <-c(10, 20, 100)
polio_Nsim <-        c(50,100, 500)

library(doParallel)
registerDoParallel()

stew(file=sprintf("pf1-%d.rda",run_level),{
  t1 <- system.time(
    pf1 <- foreach(i=1:20,.packages='pomp',
                   .options.multicore=list(set.seed=TRUE)) %dopar% try(
                     pfilter(polio,Np=polio_Np[run_level])
                   )
  )
},seed=493536993,kind="L'Ecuyer")
(L1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))

stew(sprintf("persistence-%d.rda",run_level),{
  t_sim <- system.time(
    sim <- foreach(i=1:polio_Nsim[run_level],.packages='pomp',
                   .options.multicore=list(set.seed=TRUE)) %dopar%
      simulate(polio)
  )
},seed=493536993,kind="L'Ecuyer")

no_cases_data <- sum(obs(polio)==0)
no_cases_sim <- sum(sapply(sim,obs)==0)/length(sim)
fadeout1_sim <- sum(sapply(sim,function(po)states(po)["IB",]+states(po)["IO",]<1))/length(sim)
fadeout100_sim <- sum(sapply(sim,function(po)states(po)["IB",]+states(po)["IO",]<100))/length(sim)
imports_sim <- coef(polio)["psi"]*mean(sapply(sim,function(po) mean(states(po)["SO",]+states(po)["SB1",]+states(po)["SB2",]+states(po)["SB3",]+states(po)["SB4",]+states(po)["SB5",]+states(po)["SB6",])))/12

mle_simulation <- simulate(polio,seed=127)
plot(mle_simulation)

polio_rw.sd_rp <- 0.02
polio_rw.sd_ivp <- 0.2
polio_cooling.fraction.50 <- 0.5

stew(sprintf("mif-%d.rda",run_level),{
  t2 <- system.time({
    m2 <- foreach(i=1:polio_Nreps_local[run_level],
                  .packages='pomp', .combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar% try(
                    mif2(polio,
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

    lik_m2 <- foreach(i=1:polio_Nreps_local[run_level],.packages='pomp',
                      .combine=rbind,.options.multicore=list(set.seed=TRUE)) %dopar%
                      {
                        logmeanexp(
                          replicate(polio_Nreps_eval[run_level],
                                    logLik(pfilter(polio,params=coef(m2[[i]]),Np=polio_Np[run_level]))
                          ),
                          se=TRUE)
                      }
  })
},seed=318817883,kind="L'Ecuyer")

r2 <- data.frame(logLik=lik_m2[,1],logLik_se=lik_m2[,2],t(sapply(m2,coef)))
if (run_level>1)
  write.table(r2,file="polio_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r2$logLik,digits=5)

pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(r2,logLik>max(logLik)-20))

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

stew(file=sprintf("box_eval-%d.rda",run_level),{
  t3 <- system.time({
    m3 <- foreach(i=1:polio_Nreps_global[run_level],.packages='pomp',.combine=c,
                  .options.multicore=list(set.seed=TRUE)) %dopar%
      mif2(
        m2[[1]],
        start=c(apply(polio_box,1,function(x)runif(1,x[1],x[2])),polio_fixed_params)
      )

    lik_m3 <- foreach(i=1:polio_Nreps_global[run_level],.packages='pomp',.combine=rbind,
                      .options.multicore=list(set.seed=TRUE)) %dopar% {
                        set.seed(87932+i)
                        logmeanexp(
                          replicate(polio_Nreps_eval[run_level],
                                    logLik(pfilter(polio,params=coef(m3[[i]]),Np=polio_Np[run_level]))
                          ),
                          se=TRUE)
                      }
  })
},seed=290860873,kind="L'Ecuyer")


r3 <- data.frame(logLik=lik_m3[,1],logLik_se=lik_m3[,2],t(sapply(m3,coef)))
if(run_level>1) write.table(r3,file="polio_params.csv",append=TRUE,col.names=FALSE,row.names=FALSE)
summary(r3$logLik,digits=5)

pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(r3,logLik>max(logLik)-20))

nb_lik <- function(theta) -sum(dnbinom(as.vector(obs(polio)),size=exp(theta[1]),prob=exp(theta[2]),log=TRUE))
nb_mle <- optim(c(0,-5),nb_lik)
-nb_mle$value

log_y <- log(as.vector(obs(polio))+1)
arma_fit <- arima(log_y,order=c(2,0,2),seasonal=list(order=c(1,0,1),period=12))
arma_fit$loglik-sum(log_y)

polio_params <- read.table("polio_params.csv",row.names=NULL,header=TRUE)
pairs(~logLik+psi+rho+tau+sigma_dem+sigma_env,data=subset(polio_params,logLik>max(logLik)-20))

plot(logLik~rho,data=subset(r3,logLik>max(r3$logLik)-10),log="x")

class(m3)

class(m3[[1]])

plot(m3[r3$logLik>max(r3$logLik)-10])

loglik_convergence <- do.call(cbind,conv.rec(m3[r3$logLik>max(r3$logLik)-10],"loglik"))
matplot(loglik_convergence,type="l",lty=1,ylim=max(loglik_convergence,na.rm=T)+c(-10,0))
