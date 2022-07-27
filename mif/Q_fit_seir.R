params <-
list(prefix = "Q_fit_seir")

library(tidyverse)
library(pomp)
library(doParallel)
library(doRNG)
if (.Platform$OS.type=="windows")
  options(pomp_cdir="./tmp")

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

measSIR %>%
  pomp(
    rprocess=euler(seir_step,delta.t=1/7),
    rinit=seir_rinit,
    partrans=parameter_trans(
      log=c("Beta","mu_EI"),
      logit=c("eta","rho")
    ),
    paramnames=c("N","Beta","mu_EI","mu_IR","eta","k","rho"),
    statenames=c("S","E","I","R","H")
  ) -> measSEIR

read_csv("measles_params.csv") %>%
  filter(abs(mu_IR-2)<0.001) %>%
  filter(loglik==max(loglik)) %>%
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

measSEIR %>%
  pfilter(Np=1000) -> pf1
logLik(pf1)
plot(pf1)

ncpu <- min(detectCores()-1,15)
options(cores=ncpu)
registerDoParallel()


bake(file="Q_fit_seir_local_search.rds",{
  registerDoRNG(482947940)
  foreach(i=seq_len(ncpu),.combine=c) %dopar% {
    library(tidyverse)
    library(pomp)
    measSEIR %>%
      mif2(
        Np=1000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(Beta=0.02, rho=0.02, eta=ivp(0.02), mu_EI=0.02)
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- getDoParWorkers()
  mifs_local
}) -> mifs_local

  mifs_local %>%
    traces(pars=c("loglik","Beta","mu_EI","rho","eta")) %>%
    melt() %>%
    ggplot(aes(x=iteration,y=value,group=L1,color=factor(L1)))+
    geom_line()+
    guides(color="none")+
    facet_wrap(~variable,scales="free_y")


bake(file="Q_fit_seir_lik_local.rds",{
  registerDoRNG(901242057)
  foreach(mf=mifs_local,.combine=rbind) %dopar% {
    library(tidyverse)
    library(pomp)
    evals <- replicate(10, logLik(pfilter(mf,Np=2000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2],Np=2000,nfilt=10)
  } -> local_logliks
  attr(local_logliks,"ncpu") <- getDoParWorkers()
  local_logliks
}) -> local_logliks

local_logliks$loglik
local_logliks$loglik.se

etime <- attr(mifs_local,"system.time")+attr(local_logliks,"system.time")
ncpu <- min(attr(mifs_local,"ncpu"),attr(local_logliks,"ncpu"))
mif_work <- sum(sapply(mifs_local,slot,"Nmif")*apply(sapply(mifs_local,slot,"Np"),2,mean))
pfilter_work <- with(local_logliks,sum(Np*nfilt))
efactor <- unname(etime[3]*ncpu/(mif_work+pfilter_work)*1000)
efactor

unit_cost <- (100*1000+10*2000)/1000*efactor
budget <- 60*ncpu
budget/unit_cost

freeze(
  runif_design(
    lower=c(Beta=5,rho=0.2,eta=0,mu_EI=1/3),
    upper=c(Beta=80,rho=0.9,eta=1,mu_EI=3),
    nseq=45
  ),
  seed=2062379496
)-> guesses



bake(file="Q_fit_seir_global1.rds",{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(tidyverse)
    library(pomp)
    measSEIR %>%
      mif2(
        Nmif=100, Np=1000,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(Beta=0.02, rho=0.02, eta=ivp(0.02), mu_EI=0.02),
        params=c(unlist(guess),fixed_params)
      ) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=2000) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> global1
  attr(global1,"ncpu") <- getDoParWorkers()
  global1
}) -> global1

pairs(
  ~loglik+Beta+eta+rho+mu_EI,
  data=filter(global1,loglik>max(loglik)-10),
  pch=16
)

unit_cost <- (200*1000+10*2000)/1000*efactor
budget <- 5*60*ncpu
budget/unit_cost

freeze(
  runif_design(
    lower=c(Beta=5,rho=0.2,eta=0,mu_EI=1/3),
    upper=c(Beta=80,rho=0.9,eta=1,mu_EI=3),
    nseq=150
  ),
  seed=2062379496
)-> guesses



bake(file="Q_fit_seir_global2.rds",{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(tidyverse)
    library(pomp)
    measSEIR %>%
      mif2(
        Nmif=100, Np=1000,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(Beta=0.02, rho=0.02, eta=ivp(0.02), mu_EI=0.02),
        params=c(unlist(guess),fixed_params)
      ) %>%
      continue(
        cooling.fraction=0.1
      ) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=2000) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> global2
  attr(global2,"ncpu") <- getDoParWorkers()
  global2
}) -> global2

pairs(
  ~loglik+Beta+eta+rho+mu_EI,
  data=filter(global2,loglik>max(loglik)-10),
  pch=16
)

unit_cost <- (200*1000+10*2000)/1000*efactor
budget <- 30*60*ncpu
budget/unit_cost



freeze(
  runif_design(
    lower=c(Beta=5,rho=0.2,eta=0,mu_EI=1/3),
    upper=c(Beta=80,rho=0.9,eta=1,mu_EI=3),
    nseq=1000
  ),
  seed=2062379496
)-> guesses

## registerDoRNG(1270401374)
## foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
##   library(tidyverse)
##   library(pomp)
##   measSEIR %>%
##     mif2(
##       Nmif=100, Np=1000,
##       cooling.fraction.50=0.5,
##       rw.sd=rw.sd(Beta=0.02, rho=0.02, eta=ivp(0.02), mu_EI=0.02),
##       params=c(unlist(guess),fixed_params)
##     ) %>%
##     continue(
##       cooling.fraction=0.1
##     ) -> mf
##   replicate(
##     10,
##     mf %>% pfilter(Np=2000) %>% logLik()
##   ) %>%
##     logmeanexp(se=TRUE) -> ll
##   mf %>% coef() %>% bind_rows() %>%
##     bind_cols(loglik=ll[1],loglik.se=ll[2])
## } -> global3



global3 %>%
  count(
    finite=is.finite(loglik),
    precise=loglik.se>0.2,
    high=loglik>max(loglik,na.rm=TRUE)-10
  ) %>%
  as.data.frame()

global3 %>%
  filter(
    is.finite(loglik),
    loglik>max(loglik,na.rm=TRUE)-10
  ) -> global3

pairs(
  ~loglik+Beta+eta+rho+mu_EI,
  data=global3,
  pch=16,cex=0.5
)

bind_rows(
  global1,
  global2,
  global3
) %>%
  filter(
    is.finite(loglik),
    loglik.se < 0.2
  ) -> mle_seir

read_csv("measles_params.csv") %>%
  filter(abs(mu_IR-2)<0.001) %>%
  filter(loglik==max(loglik)) -> mle_sir

