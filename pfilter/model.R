## ----model-construct----------------------------------------------------------
library(tidyverse)
library(pomp)

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;"
)

sir_init <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;"
)

dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,rho*H,give_log);"
)

rmeas <- Csnippet("
  reports = rnbinom_mu(k,rho*H);"
)

emeas <- Csnippet("
  E_reports = rho*H;"
)

read_csv("https://kingaa.github.io/sbied/pfilter/Measles_Consett_1948.csv") |>
  select(week,reports=cases) |>
  filter(week<=42) |>
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    emeasure=emeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","mu_IR","eta","rho","k","N"),
    params=c(Beta=15,mu_IR=0.5,rho=0.5,k=10,eta=0.06,N=38000)
  ) -> measSIR
