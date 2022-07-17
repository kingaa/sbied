params <-
list(prefix = "Q_fit_seir")

library(pomp)
library(tidyverse)
library(doParallel)
library(doRNG)

source("https://kingaa.github.io/sbied/pfilter/model.R")

seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

seir_rinit <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

measSEIR <- pomp(measSIR,
  rprocess=euler(seir_step,delta.t=1/7),
  rinit=seir_rinit,
  paramnames=c("N","Beta","mu_EI","mu_IR","eta","k","rho"),
  partrans=parameter_trans(
        log=c("Beta","mu_EI","mu_IR","k"),
        logit=c("eta","rho")
  ),
  statenames=c("S","E","I","R","H")
)

read_csv("measles_params.csv") %>%
  filter(
    loglik==max(loglik),
    abs(mu_IR-2)<0.001
    ) %>%
  select(-loglik,-loglik.se) -> coef(measSEIR)

coef(measSEIR,"mu_EI") <- 0.8
fixed_params <- coef(measSEIR,c("N","mu_IR","k"))
coef(measSEIR)

set.seed(1014406)
measSEIR %>%
  simulate(nsim=20,format="data.frame",include=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=(.id=="data")))+
  geom_line()+
  guides(color="none")+
  theme_bw()

pf1 <- pfilter(measSEIR,1000)
plot(pf1)
logLik(pf1)

## ncpu <- detectCores()
## options(cores=ncpu)
## registerDoParallel()

if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}

## registerDoRNG(482947940)
## foreach(i=1:20,.combine=c) %dopar% {
##   library(pomp)
##   library(tidyverse)
##   measSEIR %>%
##     mif2(
##       Np=2000, Nmif=50,
##       cooling.fraction.50=0.5,
##       rw.sd=rw.sd(Beta=0.02, rho=0.02, eta=ivp(0.02),mu_EI=0.02)
##     )
## } -> mifs_local
bake(file="Q_fit_seir_local_search.rds",{
  registerDoRNG(482947940)
  foreach(i=1:20,.combine=c) %dopar% {
    library(pomp)
    library(tidyverse)
    measSEIR %>%
      mif2(
        Np=2000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(Beta=0.02, rho=0.02, eta=ivp(0.02),mu_EI=0.02)
      )
  } -> mifs_local
}) -> mifs_local

sapply(mifs_local,logLik)

## registerDoRNG(900242057)
## foreach(mf=mifs_local,.combine=rbind) %dopar% {
##   library(pomp)
##   library(tidyverse)
##   evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
##   ll <- logmeanexp(evals,se=TRUE)
##   mf %>% coef() %>% bind_rows() %>%
##     bind_cols(loglik=ll[1],loglik.se=ll[2])
## } -> local_logliks
bake(file="Q_fit_seir_lik_local.rds",{
  registerDoRNG(900242057)
  foreach(mf=mifs_local,.combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> local_logliks
}) -> local_logliks

local_logliks$loglik

set.seed(2062379496)

runif_design(
  lower=c(Beta=5,rho=0.2,eta=0,mu_EI=1/3),
  upper=c(Beta=80,rho=0.9,eta=1,mu_EI=3),
  nseq=200
) -> guesses
mf1 <- mifs_local[[1]]

## registerDoRNG(1270401374)
## foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
##   library(pomp)
##   library(tidyverse)
##   mf1 %>%
##     mif2(params=c(unlist(guess),fixed_params),Np=1000) %>%
##     mif2(Nmif=100) -> mf
##   replicate(
##     10,
##     mf %>% pfilter(Np=2000) %>% logLik()
##   ) %>%
##     logmeanexp(se=TRUE) -> ll
##   mf %>% coef() %>% bind_rows() %>%
##     bind_cols(loglik=ll[1],loglik.se=ll[2])
## } -> results

bake(file="Q_fit_seir_global_search.rds",{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(params=c(unlist(guess),fixed_params),Np=1000) %>%
      mif2(Nmif=100) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=2000) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
}) %>%
  filter(is.finite(loglik)) -> results

pairs(~loglik+Beta+eta+rho+mu_EI,
      data=filter(results,loglik>max(loglik)-10))
