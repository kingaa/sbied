params <-
list(prefix = "Q_fit_all")

set.seed(1350254336)

library(pomp)
library(tidyverse)
library(doParallel)
library(doRNG)
registerDoParallel()

source("https://kingaa.github.io/sbied/pfilter/model.R")

read_csv("measles_params.csv") %>%
  filter(!is.na(loglik), loglik.se<1) %>%
  filter(loglik==max(loglik)) %>%
  select(-loglik,-loglik.se) -> coef(measSIR)

coef(measSIR) %>% mysignif(3) %>% t() %>% t() %>% knitr::kable()



## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.
bake(file="fitall_local_search.rds",{
  registerDoRNG(482942941)
  foreach(i=1:8,.combine=c) %dopar% {
    library(pomp)
    library(tidyverse)
    measSIR %>%
      mif2(
        Np=2000, Nmif=40,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(Beta=0.02, rho=0.02, mu_IR=0.02, k=0.02, eta=ivp(0.05)),
        partrans=parameter_trans(log=c("Beta","k","mu_IR"),logit=c("rho","eta")),
        paramnames=c("Beta","rho","k","eta","mu_IR")
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- getDoParWorkers()
  mifs_local
}) -> mifs_local
t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

mifs_local %>%
  traces() %>%
  melt() %>%
  filter(variable != "N") %>%
  ggplot(aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~variable,scales="free_y")



## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.
bake(file="fitall_lik_local.rds",{
  registerDoRNG(908222057)
  foreach(mf=mifs_local,.combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    evals <- replicate(10, logLik(pfilter(mf,Np=2000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")

pairs(~loglik+k+Beta+eta+rho+mu_IR,data=results,pch=16)

read_csv("measles_params.csv") %>%
  bind_rows(results) %>%
  arrange(-loglik) %>%
  write_csv("fitall_params.csv")

if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}

freeze(
  runif_design(
    lower=c(Beta=5,rho=0.2,eta=0,k=1,mu_IR=0.5),
    upper=c(Beta=80,rho=0.9,eta=1,k=30,mu_IR=2),
    nseq=500
  ),
  seed=2062379496
) -> guesses

mf1 <- mifs_local[[1]]
fixed_params <- coef(measSIR,c("N"))



## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.
bake(file="fitall_global_search.rds",
  dependson="guesses",{
  registerDoRNG(274481374)
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(
        Nmif=100,Np=2000,
        params=c(unlist(guess),fixed_params)
      ) %>%
      mif2(Nmif=100) %>%
      mif2(Nmif=100) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=2000) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")

read_csv("fitall_params.csv") %>%
  bind_rows(results) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv("fitall_params.csv")

results %>%
  filter(is.finite(loglik)) %>%
  filter(loglik>max(loglik)-50) %>%
  bind_rows(guesses) %>%
  mutate(type=if_else(is.na(loglik),"guess","result")) %>%
  arrange(type) -> all

pairs(~loglik+k+Beta+eta+rho+mu_IR, data=all, pch=16, cex=0.3,
      col=ifelse(all$type=="guess",grey(0.5),"red"))

read_csv("fitall_params.csv") %>%
  filter(loglik>max(loglik)-50) -> all

all %>%
  filter(loglik>max(loglik)-10) %>%
  ggplot(aes(x=k,y=loglik))+
  geom_point()+
  labs(
    x=expression(k),
    y=expression(log~L),
    title="poor man's profile likelihood"
  )

all %>%
  filter(loglik>max(loglik)-10) %>%
  ggplot(aes(x=eta,y=loglik,color=k))+
  geom_point()+
  labs(
    x=expression(eta),
    y=expression(log~L),
    title="poor man's profile likelihood"
  )

pairs(~loglik+k+Beta+eta+rho+mu_IR, pch=16, cex=0.3,
  data=filter(all,loglik>max(loglik)-10),
  col=ifelse(round(all$k,2)==10,"blue","red"))

read_csv("fitall_params.csv") %>%
  filter(loglik>max(loglik)-10,loglik.se<1) %>%
  sapply(range) -> box

freeze(
  profile_design(
    eta=seq(0.01,0.99,length=40),
    lower=box[1,c("Beta","rho","k","mu_IR")],
    upper=box[2,c("Beta","rho","k","mu_IR")],
    nprof=25, type="runif"
  ),
  seed=1893696051
)-> guesses



## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.
bake(file="fitall_eta_profile.rds",
  dependson="guesses",{
  registerDoRNG(830007657)
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(
        Nmif=100,Np=2000,
        params=c(unlist(guess),fixed_params),
        rw.sd=rw.sd(Beta=0.02, rho=0.02, mu_IR=0.02, k=0.02),
        partrans=parameter_trans(log=c("Beta","mu_IR","k"),logit="rho"),
        paramnames=c("Beta","mu_IR","k","rho")
      ) %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=2000) %>% logLik()
    ) %>% logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_eta <- attr(results,"system.time")
ncpu_eta <- attr(results,"ncpu")

results %>%
  filter(
    is.finite(loglik),
    loglik.se<1
  ) %>%
  group_by(eta) %>%
  filter(rank(-loglik)<=2) %>%
  ungroup() %>%
  mutate(
    `R[0]`=Beta/mu_IR,
    beta=Beta,
    `log~L`=loglik,
    `mu[IR]`=mu_IR
  ) %>%
  gather(variable,value,`log~L`,`mu[IR]`,rho,beta,k,`R[0]`) %>%
  ggplot(aes(x=eta,y=value))+
  geom_point()+
  labs(y=NULL,x=expression(eta))+
  facet_wrap(~variable,scales="free_y",labeller=label_parsed)

read_csv("fitall_params.csv") -> all
all %>%
  filter(loglik==max(loglik)) %>%
  pull(loglik) %>%
  round(1) -> ml
all %>%
  filter(round(k,4)==10) %>%
  filter(loglik==max(loglik)) %>%
  pull(loglik) %>%
  round(1) -> ml_constr
all %>%
  group_by(etacut=round(eta,2)) %>%
  filter(loglik==max(loglik)) %>%
  ungroup() %>%
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) %>%
  select(-loglik,-loglik.se,-etacut) %>%
  mutate(R0=Beta/mu_IR) %>%
  gather(parameter,val) %>%
  group_by(parameter) %>%
  summarize(min=min(val),max=max(val)) %>%
  mutate(min=signif(min,2),max=signif(max,2)) %>%
  gather(var,val,min,max) %>%
  spread(parameter,val) %>%
  arrange(eta) %>%
  column_to_rownames("var") %>%
  collect() -> cis
