##' Codes constructing the Martinez-Bakker et al. (2015) polio transmission model

library(pomp)

## ----data----------------------------------------------------------------
polio_data <- read.table("http://kingaa.github.io/sbied/polio/polio_wisconsin.csv")

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
