library(tidyverse)
library(pomp)
set.seed(1350254336)

source("https://kingaa.github.io/sbied/pfilter/model.R")



measSIR |>
  pfilter(Np=1000) -> pf


fixed_params <- c(N=38000, mu_IR=2, k=10)
coef(measSIR,names(fixed_params)) <- fixed_params

library(foreach)
library(doFuture)
plan(multisession)



tic <- Sys.time()
foreach(i=1:10,.combine=c,
  .options.future=list(seed=TRUE)
) %dofuture% {
  measSIR |> pfilter(Np=5000)
} -> pf

pf |> logLik() |> logmeanexp(se=TRUE) -> L_pf
L_pf
toc <- Sys.time()

pf[[1]] |> coef() |> bind_rows() |>
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) |>
  write_csv("measles_params.csv")


## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.
bake(file="local_search.rds",{
  foreach(i=1:20,.combine=c,
    .options.future=list(seed=482947940)
  ) %dofuture% {
    measSIR |>
      mif2(
        Np=2000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, rho=0.02, eta=ivp(0.02)),
        partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
        paramnames=c("Beta","rho","eta")
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- nbrOfWorkers()
  mifs_local
}) -> mifs_local
t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")

mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")



bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind,
    .options.future=list(seed=900242057)
  ) %dofuture% {
    evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- nbrOfWorkers()
  results
}) -> results
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")

pairs(~loglik+Beta+eta+rho,data=results,pch=16)

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  arrange(-loglik) |>
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


bake(file="global_search.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=1270401374)
    ) %dofuture% {
      mf1 |>
        mif2(params=c(guess,fixed_params)) |>
        mif2(Nmif=100) -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()
      ) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results
t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-50) |>
  bind_rows(guesses) |>
  mutate(type=if_else(is.na(loglik),"guess","result")) |>
  arrange(type) -> all

pairs(~loglik+Beta+eta+rho, data=all, pch=16, cex=0.3,
  col=ifelse(all$type=="guess",grey(0.5),"red"))

all |>
  filter(type=="result") |>
  filter(loglik>max(loglik)-10) |>
  ggplot(aes(x=eta,y=loglik))+
  geom_point()+
  labs(
    x=expression(eta),
    title="poor man's profile likelihood"
  )

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box
box

freeze(seed=1196696958,
  profile_design(
    eta=seq(0.01,0.95,length=40),
    lower=box[1,c("Beta","rho")],
    upper=box[2,c("Beta","rho")],
    nprof=15, type="runif"
  )) -> guesses
plot(guesses)



mf1 <- mifs_local[[1]]
bake(file="eta_profile.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=830007657)
    ) %dofuture% {
      mf1 |>
        mif2(params=c(guess,fixed_params),
          rw.sd=rw_sd(Beta=0.02,rho=0.02)) |>
        mif2(Nmif=100,cooling.fraction.50=0.3) -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) -> results
t_eta <- attr(results,"system.time")
ncpu_eta <- attr(results,"ncpu")

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-10) -> all

pairs(~loglik+Beta+eta+rho,data=all,pch=16)

results |>
  ggplot(aes(x=eta,y=loglik))+
  geom_point()

results |>
  filter(is.finite(loglik)) |>
  group_by(round(eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  filter(loglik>max(loglik)-20) |>
  ggplot(aes(x=eta,y=loglik))+
  geom_point()

maxloglik <- max(results$loglik,na.rm=TRUE)
ci.cutoff <- maxloglik-0.5*qchisq(df=1,p=0.95)

results |>
  filter(is.finite(loglik)) |>
  group_by(round(eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  ggplot(aes(x=eta,y=loglik))+
  geom_point()+
  geom_smooth(method="loess",span=0.25)+
  geom_hline(color="red",yintercept=ci.cutoff)+
  lims(y=maxloglik-c(5,0))

results |>
  filter(is.finite(loglik)) |>
  group_by(round(eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  mutate(in_ci=loglik>max(loglik)-1.92) |>
  ggplot(aes(x=eta,y=rho,color=in_ci))+
  geom_point()+
  labs(
    color="inside 95% CI?",
    x=expression(eta),
    y=expression(rho),
    title="profile trace"
  )

results |>
  filter(is.finite(loglik)) |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(rho),max=max(rho)) -> rho_ci

read_csv("measles_params.csv") |>
  group_by(cut=round(rho,2)) |>
  filter(rank(-loglik)<=10) |>
  ungroup() |>
  arrange(-loglik) |>
  select(-cut,-loglik,-loglik.se) -> guesses


mf1 <- mifs_local[[1]]
bake(file="rho_profile.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=2105684752)
    ) %dofuture% {
      mf1 |>
        mif2(params=guess,
          rw.sd=rw_sd(Beta=0.02,eta=ivp(0.02))) |>
        mif2(Nmif=100,cooling.fraction.50=0.3) |>
        mif2() -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) -> results
t_rho <- attr(results,"system.time")
ncpu_rho <- attr(results,"ncpu")
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

results |>
  filter(is.finite(loglik)) -> results

pairs(~loglik+Beta+eta+rho,data=results,pch=16)

results |>
  filter(loglik>max(loglik)-10,loglik.se<1) |>
  group_by(round(rho,2)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  ggplot(aes(x=rho,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )

results |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(rho),max=max(rho)) -> rho_ci

freeze(seed=55266255,
  runif_design(
    lower=c(Beta=5,mu_IR=0.2,eta=0),
    upper=c(Beta=80,mu_IR=5,eta=0.99),
    nseq=1000
  )) |>
  mutate(
    rho=0.6,
    k=10,
    N=38000
  ) -> guesses





bake(file="global_search2.rds",
  dependson=guesses,{
    foreach(
      guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=610408798)
    ) %dofuture% {
      measSIR |>
        mif2(params=guess, Np=2000, Nmif=100,
          cooling.fraction.50=0.5,
          partrans=parameter_trans(
            log=c("Beta","mu_IR"),
            logit="eta"), paramnames=c("Beta","mu_IR","eta"),
          rw.sd=rw_sd(Beta=0.02,mu_IR=0.02,eta=ivp(0.02))) -> mf
      mf |>
        mif2(
          Nmif=100,rw.sd=rw_sd(Beta=0.01,mu_IR=0.01,eta=ivp(0.01))
        ) |>
        mif2(
          Nmif=100,
          rw.sd=rw_sd(Beta=0.005,mu_IR=0.005,eta=ivp(0.005))
        ) -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()
      ) |> logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results
t_expglob <- attr(results,"system.time")
ncpu_expglob <- attr(results,"ncpu")

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20) -> all

pairs(~loglik+rho+mu_IR+Beta+eta,data=all,pch=16,cex=0.3,
  col=if_else(round(all$rho,3)==0.6,1,4))

results |>
  filter(loglik>max(loglik)-20,loglik.se<1) |>
  ggplot(aes(x=mu_IR,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )

read_csv("measles_params.csv") |>
  filter(
    loglik>max(loglik)-20,
    loglik.se<2,
    abs(rho-0.6)<0.01
  ) |>
  sapply(range) -> box

freeze(seed=610408798,
  profile_design(
    mu_IR=seq(0.2,2,by=0.1),
    lower=box[1,c("Beta","eta")],
    upper=box[2,c("Beta","eta")],
    nprof=100, type="runif"
  )) |>
  mutate(
    N=38000,
    rho=0.6,
    k=10
  ) -> guesses



bake(file="mu_IR_profile1.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=610408798)
    ) %dofuture% {
      measSIR |>
        mif2(params=guess, Np=2000, Nmif=100,
          partrans=parameter_trans(log="Beta",logit="eta"),
          paramnames=c("Beta","eta"), cooling.fraction.50=0.5,
          rw.sd=rw_sd(Beta=0.02,eta=ivp(0.02))
        ) |>
        mif2(Nmif=100) |>
        mif2(Nmif=100,rw.sd=rw_sd(Beta=0.01,eta=ivp(0.01))) |>
        mif2(Nmif=100,rw.sd=rw_sd(Beta=0.005,eta=ivp(0.005))) -> mf
      replicate(10,mf |> pfilter(Np=5000) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results
t_muIR <- attr(results,"system.time")
ncpu_muIR <- attr(results,"ncpu")

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

results |>
  group_by(round(mu_IR,2)) |>
  filter(rank(-loglik)<=1) |>
  ungroup() |>
  ggplot(aes(x=mu_IR,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )
