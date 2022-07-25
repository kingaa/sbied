params <-
list(prefix = "windows")

set.seed(1350254336)

dir.create("tmp")
options(pomp_cdir="./tmp")

library(tidyverse)
library(pomp)

library(doParallel)
d <- detectCores()
registerDoParallel(d-1)

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
    times="week", t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","mu_IR","eta","rho","k","N"),
    params=c(Beta=15,mu_IR=0.5,rho=0.5,k=10,eta=0.06,N=38000)
  ) -> measSIR


foreach (i=1:4,.combine=c) %dopar% {
  library(pomp)
  measSIR %>%
    mif2(
      Np=1000, Nmif=5,
      cooling.fraction.50=0.5,
      rw.sd=rw.sd(Beta=0.2, rho=0.2, eta=ivp(0.2)),
      partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
      paramnames=c("Beta","rho","eta")
    )
} -> mifs_local

## Compilation is done here, since the call to `parameter_trans`
## leads to compilation of a C snippet.
measSIR %>%
  pomp(
    partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
    paramnames=c("Beta","rho","eta")
  ) -> measSIR2

## No compilation is triggered inside the parallel code block.
foreach (i=1:4,.combine=c) %dopar% {
  library(pomp)
  measSIR2 %>%
    mif2(
      Np=1000, Nmif=5,
      cooling.fraction.50=0.5,
      rw.sd=rw.sd(Beta=0.2, rho=0.2, eta=ivp(0.2))
    )
} -> mifs_local
