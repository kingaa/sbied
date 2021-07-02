params <-
list(prefix = "pfnm/")

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
  pomp(partrans=parameter_trans(log=c("Beta","mu_IR"),logit=c("rho","eta")),
    paramnames=c("Beta","mu_IR","eta","rho")) -> measSIR

coef(measSIR)

neg.ll <- function (par, est) {
  try(
    freeze({
      allpars <- coef(measSIR,transform=TRUE)
      allpars[est] <- par
      theta <- partrans(measSIR,allpars,dir="fromEst")
      pfilter(measSIR,params=theta,Np=2000)
    },
    seed=915909831
    )
  ) -> pf
  if (inherits(pf,"try-error")) 1e10 else -logLik(pf)
}

## use Nelder-Mead with fixed RNG seed
estpars <- c("Beta","mu_IR","eta")
optim(
  par=coef(measSIR,estpars,transform=TRUE),
  est=estpars,
  fn=neg.ll,
  method="Nelder-Mead",
  control=list(maxit=400,trace=0)
) -> fit

mle <- measSIR
coef(mle,estpars,transform=TRUE) <- fit$par
coef(mle)

fit$val

lls <- replicate(n=5,logLik(pfilter(mle,Np=20000)))
ll <- logmeanexp(lls,se=TRUE); ll

mle %>% simulate(nsim=10,format="data.frame",include.data=TRUE) -> sims

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  guides(color=FALSE)+
  geom_line()
