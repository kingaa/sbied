## ----opts,include=FALSE--------------------------------------------------
library(pomp)
library(knitr)
prefix <- "contacts"
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
#  out.width="\\linewidth",
  dpi=300,
  dev='png',
  dev.args=list(bg='transparent')
  )

options(
  keep.source=TRUE,
  encoding="UTF-8"
  )

library(ggplot2)
theme_set(theme_bw())

## ----data----------------------------------------------------------------
contact_data <- read.table(file="contacts.csv",header=TRUE)
matplot(t(contact_data[1:15,1:4]),
  ylab="total sexual contacts",xlab="6-month intervals", 
  type="l",xaxp=c(1,4,3))

## ----package-------------------------------------------------------------
require(pomp)
packageVersion("pomp")

## ----statenames----------------------------------------------------------
contact_statenames = c("X","D","R","C1","C2","C3","C4")
contact_obsnames = c("y1", "y2", "y3", "y4")

## ----params--------------------------------------------------------------
contact_paramnames = c("mu_X","sigma_X","mu_D","sigma_D","mu_R","sigma_R","alpha")
contact_mle <- c(
    mu_X = 1.75,
    sigma_X = 2.67,
    mu_D = 3.81,
    sigma_D = 4.42,
    mu_R = 0.04,
    sigma_R = 0,
    alpha = 0.90
  )

## ----rprocess------------------------------------------------------------
contact_rprocess <- Csnippet("
  double C[4] = {0,0,0,0}, tol=0.000001; 
  int j;
  double Z, Zcum;

  D = (sigma_D < tol || mu_D < tol) ? mu_D : rgamma( pow(mu_D, 2) / pow(sigma_D, 2) , pow(sigma_D, 2) / mu_D );
  R = (sigma_R < tol || mu_R < tol) ? mu_R : rgamma( pow(mu_R, 2) / pow(sigma_R, 2) , pow(sigma_R, 2) / mu_R );
  X = (sigma_X < tol || mu_X < tol) ? mu_X : rgamma( pow(mu_X, 2) / pow(sigma_X, 2) , pow(sigma_X, 2) / mu_X );
  Z = (R < tol) ? 1/tol : rexp(1/R);
  for(j=0;j<4;j++){
    Zcum = Z;
    while(Zcum < 6){
       C[j] += Z * X;
       Z = (R < tol) ? 1/tol : rexp(1/R);
       X = (sigma_X < tol || mu_X < tol) ? mu_X : rgamma( pow(mu_X, 2) / pow(sigma_X, 2) , pow(sigma_X, 2) / mu_X );
       Zcum += Z;
    }
    C[j] += (6 - (Zcum - Z)) * X;
    C[j] *= pow(alpha,j); 
    Z = Zcum - 6;
  }
  C1 = C[0]; C2 = C[1]; C3 = C[2]; C4 = C[3];
")

## ----measure-------------------------------------------------------------
contact_dmeasure <- Csnippet("
  double y[4] = {y1,y2,y3,y4};
  double C[4] = {C1,C2,C3,C4};
  int j;
  lik = 0;
  for(j=0;j<4;j++) lik += dnbinom(y[j], D, D/(D+C[j]), 1);
  lik = give_log ? lik : exp(lik);
")

contact_rmeasure <- Csnippet("
  double y[4];
  double C[4] = {C1,C2,C3,C4};
  int j;
  for(j=0;j<4;j++) y[j] = rnbinom(D, D/(D+C[j]));
  y1=y[0]; y2=y[1]; y3=y[2]; y4=y[3];
")

## ----pomp----------------------------------------------------------------
contacts <- pomp(
  data=contact_data,
  times="individual",
  t0=0,
  params=contact_mle,
  rprocess = discrete.time.sim(step.fun=contact_rprocess,delta.t=1),
  dmeasure = contact_dmeasure,
  rmeasure = contact_rmeasure,
  obsnames = contact_obsnames,
  statenames = contact_statenames,
  paramnames = contact_paramnames,
  initializer = function (params, t0,  ...) c(X=0,D=0,R=0,C1=0,C2=0,C3=0,C4=0)
) 

plot(contacts)

## ----simulate------------------------------------------------------------
s1 <- simulate(contacts,seed=1)
apply(obs(s1),1,mean)
apply(states(s1),1,mean)
plot(s1)

## ----pfilter-------------------------------------------------------------
require(doParallel)
cores <- 10
registerDoParallel(cores)
t1 <- system.time(
 pf1 <- foreach(i=1:cores,.packages='pomp') %dopar% try(pfilter(contacts,Np=2000,seed=2015+i))
)
(loglik1 <- sapply(pf1,logLik))

## ----transformations-----------------------------------------------------
contact_toEstimationScale <- Csnippet("
  Tmu_X = log(mu_X);
  Tsigma_X = log(sigma_X);
  Tmu_D = log(mu_D);
  Tsigma_D = log(sigma_D);
  Tmu_R = log(mu_R);
  Talpha = log(alpha/(1-alpha));
")

contact_fromEstimationScale <- Csnippet("
  Tmu_X = exp(mu_X);
  Tsigma_X = exp(sigma_X);
  Tmu_D = exp(mu_D);
  Tsigma_D = exp(sigma_D);
  Tmu_R = exp(mu_R);
  Talpha = exp(alpha)/(1+exp(alpha));
")

contacts_with_trans <- pomp(contacts,
  paramnames = contact_paramnames,
  fromEstimationScale=contact_fromEstimationScale,
  toEstimationScale=contact_toEstimationScale
)

## ----statenames2---------------------------------------------------------
contact2_statenames = c("X","D","R","C","Z")
contact2_obsnames = "y"

## ----rprocess2-----------------------------------------------------------
contact2_rprocess <- Csnippet("
  double Zcum, tol=0.000001;
  if( (int)t % 4 == 0) { 
    D = (sigma_D < tol || mu_D < tol) ? mu_D : rgamma( pow(mu_D, 2) / pow(sigma_D, 2) , pow(sigma_D, 2) / mu_D );
    R = (sigma_R < tol || mu_R < tol) ? mu_R : rgamma( pow(mu_R, 2) / pow(sigma_R, 2) , pow(sigma_R, 2) / mu_R );
    X = (sigma_X < tol || mu_X < tol) ? mu_X : rgamma( pow(mu_X, 2) / pow(sigma_X, 2) , pow(sigma_X, 2) / mu_X );
    Z = (R < tol) ? 1/tol : rexp(1/R);
  }
  C = 0;
  Zcum = Z;
  while(Zcum < 6){
       C  += Z * X;
       Z = (R < tol) ? 1/tol : rexp(1/R);
       X = (sigma_X < tol || mu_X < tol) ? mu_X : rgamma( pow(mu_X, 2) / pow(sigma_X, 2) , pow(sigma_X, 2) / mu_X );
       Zcum += Z;
  }
  C += (6 - (Zcum - Z)) * X;
  C *= pow(alpha, (int)t % 4 ); 
  Z = Zcum - 6;
")

## ----measurement2--------------------------------------------------------

contact2_dmeasure <- Csnippet("
  lik = dnbinom(y, D, D/(D+C), give_log);
")

contact2_rmeasure <- Csnippet("
  y = rnbinom(D, D/(D+C));
")

## ----pomp2---------------------------------------------------------------
contacts2 <- pomp(
  data=read.table(file="contacts2.csv",header=TRUE),
  times="obs",
  t0=0,
  params=contact_mle,
  rprocess = discrete.time.sim(step.fun=contact2_rprocess,delta.t=1),
  dmeasure = contact2_dmeasure,
  rmeasure = contact2_rmeasure,
  obsnames = contact2_obsnames,
  statenames = contact2_statenames,
  paramnames = contact_paramnames,
  fromEstimationScale=contact_fromEstimationScale,
  toEstimationScale=contact_toEstimationScale,
  initializer = function (params, t0,  ...) c(X=0,D=0,R=0,C=0,Z=0)
) 

plot(contacts2)

## ----simulateVersion2----------------------------------------------------
s2 <- simulate(contacts2,seed=2)
apply(obs(s2),1,mean)
apply(states(s2),1,mean)
plot(s2)

## ----pfilter2------------------------------------------------------------
t2 <- system.time(
 pf2 <- foreach(i=1:cores,.packages='pomp') %dopar% try(pfilter(contacts2,Np=2000,seed=2015+i))
)
(loglik2 <- sapply(pf2,logLik))

## ----mif-----------------------------------------------------------------
t3 <- system.time(
  m2 <- foreach(i=1:cores,.packages='pomp') %dopar% try( 
    mif2(contacts2,
      seed=1903+i,
      Nmif=10,
      Np=200,
      cooling.fraction.50=0.5,
      cooling.type="geometric",
      transform=TRUE,
      rw.sd=rw.sd(mu_X=0.02,
              sigma_X=0.02,
              mu_D = 0.02,
              sigma_D=0.02,
              mu_R=0.02,
              sigma_R =0.02,
              alpha=0.02)
    )
  )
)
params_new <- coef( m2[[which.max( sapply(m2,logLik) )]] )
pf3 <- foreach(i=1:cores,.packages='pomp') %dopar% try(
  pfilter(contacts2,params=params_new,Np=1000,seed=1809+i)
)
(loglik_new <- logmeanexp(sapply(pf3,logLik),se=TRUE))

