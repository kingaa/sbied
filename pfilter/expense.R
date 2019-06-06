library(tidyverse)
library(pomp)
stopifnot(packageVersion("pomp")>="2.1")
options(stringsAsFactors=FALSE)

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
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
  lik = dbinom(reports,H,rho,give_log);
")

rmeas <- Csnippet("
  reports = rbinom(H,rho);
")

read_csv("https://kingaa.github.io/sbied/pfilter/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) %>%
  filter(week<=42) %>%
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/6),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","gamma","eta","rho","N"),
    params=c(Beta=15,gamma=0.5,rho=0.5,eta=0.06,N=38000)
  ) -> measSIR

Nps <- ceiling(10^seq(1,5,by=0.2))
times <- c()
for (np in Nps) {
  times <- c(
    times,
    system.time(measSIR %>% pfilter(Np=np))[3]
  )
}

plot(Nps,times)
lm(times~Nps) -> fit
abline(fit)
summary(fit)
