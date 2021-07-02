params <-
list(prefix = "direct/")

library(plyr)
library(tidyverse)
library(pomp)
stopifnot(packageVersion("pomp")>="3.0")
set.seed(594709947L)

knitr::read_chunk("model.R")
library(tidyverse)
library(pomp)

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,rho*H,give_log);
  ")

rmeas <- Csnippet("
  reports = rnbinom_mu(k,rho*H);
  ")

read_csv("https://kingaa.github.io/sbied/pfilter/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) %>%
  filter(week<=42) %>%
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","mu_IR","eta","rho","k","N"),
    params=c(Beta=15,mu_IR=0.5,rho=0.5,k=10,eta=0.06,N=38000)
  ) -> measSIR

measSIR %>%
  simulate(nsim=5000,format="arrays") -> x
sims <- coef(measSIR,"rho")*x$states["H",,]
matplot(time(measSIR),t(sims[1:50,]),type='l',lty=1,
  xlab="time",ylab=expression(rho*H),bty='l',col='blue')
lines(time(measSIR),obs(measSIR,"reports"),lwd=2,col='black')

ell <- dmeasure(measSIR,y=obs(measSIR),x=x$states,times=time(measSIR),log=TRUE,
  params=coef(measSIR))
dim(ell)

ell <- apply(ell,1,sum)
summary(ell)
summary(exp(ell))
