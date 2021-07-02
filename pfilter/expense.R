params <-
list(prefix = "expense/")

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
stopifnot(packageVersion("pomp")>="3.1")

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

expand_grid(
  Np=ceiling(10^seq(1,5,by=0.2))
) -> design

registerDoParallel()
registerDoRNG()

foreach (
  expt=iter(design,"row"),
  .combine=bind_rows
) %dopar% {
  system.time(measSIR %>% pfilter(Np=expt$Np))[3] -> expt$time
  expt
} -> resultA

resultA %>%
  ggplot(aes(x=Np,y=time))+
  geom_point()+
  geom_smooth(method="lm")+
  expand_limits(x=0,y=0)

lm(time~Np,data=resultA) -> fit
summary(fit)

measSIR %>%
  window(end=21) -> shortMeasSIR

foreach (
  expt=iter(design,"row"),
  .combine=bind_rows
) %dopar% {
  system.time(shortMeasSIR %>% pfilter(Np=expt$Np))[3] -> expt$time
  expt
} -> resultB

bind_rows(
  long=resultA,
  short=resultB,
  .id="length"
) %>%
  mutate(
    n=case_when(
      length=="short"~length(time(shortMeasSIR)),
      length=="long"~length(time(measSIR))
    )
  ) -> result

result %>%
  ggplot(aes(x=time,y=Np,group=n,color=factor(n)))+
  geom_point()+
  labs(color="n")+
  geom_smooth(method="lm")+
  expand_limits(x=0,y=0)

lm(time~n*Np,data=result) -> fit
summary(fit)
