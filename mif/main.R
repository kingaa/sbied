library(tidyverse)
library(pomp)
stopifnot(getRversion() >= "4.0")
stopifnot(packageVersion("pomp")>="3.0")
set.seed(1350254336)

library(tidyverse)

courseurl <- "https://kingaa.github.io/sbied/"
datafile <- "mif/Measles_Consett_1948.csv"

read_csv(paste0(courseurl,datafile)) %>%
  select(week,reports=cases) %>%
  filter(week<=42) -> dat

dat %>%
  ggplot(aes(x=week,y=reports))+
  geom_line()



library(pomp)

sir_step <- Csnippet("
double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
S -= dN_SI;
I += dN_SI - dN_IR;
H += dN_IR;
")

sir_init <- Csnippet("
S = nearbyint(eta*N);
I = 1;
H = 0;
")

dmeas <- Csnippet("
lik = dnbinom_mu(reports,k,rho*H,give_log);
")

rmeas <- Csnippet("
reports = rnbinom_mu(k,rho*H);
")

dat %>%
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    partrans=parameter_trans(
      log=c("Beta"),
      logit=c("rho","eta")
    ),
    statenames=c("S","I","H"),
    paramnames=c("Beta","mu_IR","eta","rho","k","N")
  ) -> measSIR

params <- c(Beta=20,mu_IR=2,rho=0.5,k=10,eta=0.1,N=38000)

measSIR %>%
  simulate(params=params,nsim=10,format="data.frame") -> y


measSIR %>%
  pfilter(Np=1000,params=params) -> pf


fixed_params <- c(N=38000, mu_IR=2, k=10)

library(foreach)
library(doParallel)
registerDoParallel()





library(doRNG)
registerDoRNG(625904618)
tic <- Sys.time()
foreach(i=1:10,.combine=c) %dopar% {
  library(pomp)
  measSIR %>% pfilter(params=params,Np=10000)
} -> pf

pf %>% logLik() %>% logmeanexp(se=TRUE) -> L_pf
L_pf
toc <- Sys.time()

pf[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) %>%
  write_csv("measles_params.csv")


registerDoRNG(482947940)
bake(file="local_search.rds",{
  foreach(i=1:20,.combine=c) %dopar% {
    library(pomp)
    library(tidyverse)
    measSIR %>%
      mif2(
        params=params,
        Np=2000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw.sd(Beta=0.02, rho=0.02, eta=ivp(0.02))
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- getDoParWorkers()
  mifs_local
}) -> mifs_local
t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

## dat %>%
##   pomp(
##     times="week",t0=0,
##     rprocess=euler(sir_step,delta.t=1/7),
##     rinit=sir_init, rmeasure=rmeas, dmeasure=dmeas,
##     partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
##     accumvars="H", statenames=c("S","I","H"),
##     paramnames=c("Beta","mu_IR","eta","rho","N"),
##     cdir=".", cfile="measSIR"
##   ) -> measSIR

mifs_local %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~variable,scales="free_y")



registerDoRNG(900242057)
bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    evals <- replicate(10, logLik(pfilter(mf,Np=10000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")

pairs(~loglik+Beta+eta+rho,data=results,pch=16)

read_csv("measles_params.csv") %>%
  bind_rows(results) %>%
  arrange(-loglik) %>%
  write_csv("measles_params.csv")

if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}

set.seed(2062379496)

runif_design(
  lower=c(Beta=5,rho=0.2,eta=0),
  upper=c(Beta=80,rho=0.9,eta=1),
  nseq=400
) -> guesses

mf1 <- mifs_local[[1]]


bake(file="global_search.rds",{
  registerDoRNG(1270401374)
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(params=c(unlist(guess),fixed_params)) %>%
      mif2(Nmif=100) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=10000) %>% logLik()
    ) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) %>%
  filter(is.finite(loglik)) -> results
t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")
read_csv("measles_params.csv") %>%
  bind_rows(results) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv("measles_params.csv")

read_csv("measles_params.csv") %>%
  filter(loglik>max(loglik)-50) %>%
  bind_rows(guesses) %>%
  mutate(type=if_else(is.na(loglik),"guess","result")) %>%
  arrange(type) -> all

pairs(~loglik+Beta+eta+rho, data=all,
      col=ifelse(all$type=="guess",grey(0.5),"red"),pch=16)

all %>%
  filter(type=="result") %>%
  filter(loglik>max(loglik)-10) %>%
  ggplot(aes(x=eta,y=loglik))+
  geom_point()+
  labs(
    x=expression("eta"),
    title="poor man's profile likelihood"
  )

read_csv("measles_params.csv") %>%
  filter(loglik>max(loglik)-20,loglik.se<2) %>%
  sapply(range) -> box
box

set.seed(1196696958)
profile_design(
  eta=seq(0.01,0.95,length=40),
  lower=box[1,c("Beta","rho")],
  upper=box[2,c("Beta","rho")],
  nprof=15, type="runif"
) -> guesses
plot(guesses)



mf1 <- mifs_local[[1]]
registerDoRNG(830007657)
bake(file="eta_profile.rds",{
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(params=c(unlist(guess),fixed_params),
           rw.sd=rw.sd(Beta=0.02,rho=0.02)) %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=10000) %>% logLik()) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_eta <- attr(results,"system.time")
ncpu_eta <- attr(results,"ncpu")

read_csv("measles_params.csv") %>%
  bind_rows(results) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv("measles_params.csv")

read_csv("measles_params.csv") %>%
  filter(loglik>max(loglik)-10) -> all

pairs(~loglik+Beta+eta+rho,data=all,pch=16)

results %>%
  ggplot(aes(x=eta,y=loglik))+
  geom_point()

results %>%
  filter(is.finite(loglik)) %>%
  group_by(round(eta,5)) %>%
  filter(rank(-loglik)<3) %>%
  ungroup() %>%
  filter(loglik>max(loglik)-20) %>%
  ggplot(aes(x=eta,y=loglik))+
  geom_point()

maxloglik <- max(results$loglik,na.rm=TRUE)
ci.cutoff <- maxloglik-0.5*qchisq(df=1,p=0.95)

results %>%
  filter(is.finite(loglik)) %>%
  group_by(round(eta,5)) %>%
  filter(rank(-loglik)<3) %>%
  ungroup() %>%
  ggplot(aes(x=eta,y=loglik))+
  geom_point()+
  geom_smooth(method="loess",span=0.25)+
  geom_hline(color="red",yintercept=ci.cutoff)+
  lims(y=maxloglik-c(5,0))

results %>%
  filter(is.finite(loglik)) %>%
  group_by(round(eta,5)) %>%
  filter(rank(-loglik)<3) %>%
  ungroup() %>%
  mutate(in_ci=loglik>max(loglik)-1.92) %>%
  ggplot(aes(x=eta,y=rho,color=in_ci))+
  geom_point()+
  labs(
    color="inside 95% CI?",
    x=expression(eta),
    y=expression(rho),
    title="profile trace"
  )

results %>%
  filter(is.finite(loglik)) %>%
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) %>%
  summarize(min=min(rho),max=max(rho)) -> rho_ci

read_csv("measles_params.csv") %>%
  group_by(cut=round(rho,2)) %>%
  filter(rank(-loglik)<=10) %>%
  ungroup() %>%
  select(-cut,-loglik,-loglik.se) -> guesses


mf1 <- mifs_local[[1]]
registerDoRNG(2105684752)
bake(file="rho_profile.rds",{
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    mf1 %>%
      mif2(params=guess,
           rw.sd=rw.sd(Beta=0.02,eta=ivp(0.02))) %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) %>%
      mif2() -> mf
    replicate(
      10,
      mf %>% pfilter(Np=10000) %>% logLik()) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) -> results
t_rho <- attr(results,"system.time")
ncpu_rho <- attr(results,"ncpu")
read_csv("measles_params.csv") %>%
  bind_rows(results) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv("measles_params.csv")

results %>%
  filter(is.finite(loglik)) -> results

pairs(~loglik+Beta+eta+rho,data=results,pch=16)

results %>%
  filter(loglik>max(loglik)-10,loglik.se<1) %>%
  group_by(round(rho,2)) %>%
  filter(rank(-loglik)<3) %>%
  ungroup() %>%
  ggplot(aes(x=rho,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )

results %>%
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) %>%
  summarize(min=min(rho),max=max(rho)) -> rho_ci

set.seed(55266255)
runif_design(
  lower=c(Beta=5,mu_IR=0.2,eta=0),
  upper=c(Beta=80,mu_IR=5,eta=0.95),
  nseq=1000
) %>%
  mutate(
    rho=0.6,
    k=10,
    N=38000
  ) -> guesses







registerDoRNG(610408798)
bake(file="global_search2.rds",{
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    measSIR %>%
      mif2(params=guess, Np=2000, Nmif=100,
           cooling.fraction.50=0.5,
           partrans=parameter_trans(
             log=c("Beta","mu_IR"),
             logit="eta"), paramnames=c("Beta","mu_IR","eta"),
           rw.sd=rw.sd(Beta=0.02,mu_IR=0.02,eta=ivp(0.02))) -> mf
    mf %>%
      mif2(Nmif=100) %>%
      mif2(Nmif=100) %>%
      mif2(Nmif=100) -> mf
    mf %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) %>%
      mif2(Nmif=100,cooling.fraction.50=0.1) %>%
      mif2(Nmif=100,cooling.fraction.50=0.1) %>%
      mif2(Nmif=100,cooling.fraction.50=0.1) -> mf
    replicate(
      10,
      mf %>% pfilter(Np=10000) %>% logLik()
    ) %>% logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) %>%
  filter(is.finite(loglik)) -> results
t_expglob <- attr(results,"system.time")
ncpu_expglob <- attr(results,"ncpu")

read_csv("measles_params.csv") %>%
  bind_rows(results) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv("measles_params.csv")

read_csv("measles_params.csv") %>%
  filter(loglik>max(loglik)-20) -> all

pairs(~loglik+rho+mu_IR+Beta+eta,data=all,pch=16)

pairs(~loglik+rho+mu_IR+Beta+eta,data=results,pch=16)

results %>%
  filter(loglik>max(loglik)-20,loglik.se<1) %>%
  ggplot(aes(x=mu_IR,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )

read_csv("measles_params.csv") %>%
  filter(
    loglik>max(loglik)-20,
    loglik.se<2,
    abs(rho-0.6)<0.01
  ) %>%
  sapply(range) -> box

set.seed(610408798)
profile_design(
  mu_IR=seq(0.5,2,by=0.1),
  lower=box[1,c("Beta","eta")],
  upper=box[2,c("Beta","eta")],
  nprof=100, type="runif"
) %>%
  mutate(
    N=38000,
    rho=0.6,
    k=10
  ) -> guesses



registerDoRNG(610408798)
bake(file="mu_IR_profile1.rds",{
  foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
    library(pomp)
    library(tidyverse)
    measSIR %>%
      mif2(params=guess, Np=2000, Nmif=100,
           partrans=parameter_trans(log="Beta",logit="eta"),
           paramnames=c("Beta","eta"), cooling.fraction.50=0.5,
           rw.sd=rw.sd(Beta=0.02,eta=ivp(0.02))
           ) %>% mif2(Nmif=100) %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) %>%
      mif2(Nmif=100,cooling.fraction.50=0.1) -> mf
    replicate(10,mf %>% pfilter(Np=100000) %>% logLik()) %>%
      logmeanexp(se=TRUE) -> ll
    mf %>% coef() %>% bind_rows() %>%
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- getDoParWorkers()
  results
}) %>%
  filter(is.finite(loglik)) -> results
t_muIR <- attr(results,"system.time")
ncpu_muIR <- attr(results,"ncpu")

read_csv("measles_params.csv") %>%
  bind_rows(results) %>%
  filter(is.finite(loglik)) %>%
  arrange(-loglik) %>%
  write_csv("measles_params.csv")

results %>%
  group_by(round(mu_IR,2)) %>%
  filter(rank(-loglik)<3) %>%
  ungroup() %>%
  ggplot(aes(x=mu_IR,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )
