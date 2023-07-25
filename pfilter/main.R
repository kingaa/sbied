library(tidyverse)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
options(stringsAsFactors=FALSE)
set.seed(1350254336)







source("https://kingaa.github.io/sbied/pfilter/model.R")

library(pomp)
measSIR |>
  pfilter(Np=5000) -> pf
logLik(pf)

plan(multisession)
foreach (
  i=1:10,
  .combine=c,
  .options.future=list(seed=652643293)
) %dofuture% {
  measSIR |> pfilter(Np=5000)
} -> pf
logLik(pf) -> ll
logmeanexp(ll,se=TRUE)





## What is this 'bake' function?
## See https://kingaa.github.io/sbied/misc/bake.html
## for an explanation.
bake(file="like-slice.rds",{
  slice_design(
    center=coef(measSIR),
    Beta=rep(seq(from=5,to=30,length=40),each=3),
    mu_IR=rep(seq(from=0.2,to=2,length=40),each=3)
  ) -> p
  library(iterators)
  plan(multisession)
  foreach (
    theta=iter(p,"row"),
    .combine=rbind,
    .options.future=list(seed=108028909)
  ) %dofuture% {
    measSIR |> pfilter(params=theta,Np=5000) -> pf
    theta$loglik <- logLik(pf)
    theta
  } -> p
}) -> p

p |>
  pivot_longer(c(Beta,mu_IR)) |>
  filter(name==slice) |>
  ggplot(aes(x=value,y=loglik,color=name))+
  geom_point()+
  facet_grid(~name,scales="free_x")+
  guides(color="none")+
  labs(x="parameter value",color="")





bake(file="pfilter-grid1.rds",{
  expand.grid(
    Beta=rep(seq(from=10,to=30,length=40),each=3),
    mu_IR=rep(seq(from=0.4,to=1.5,length=40),each=3),
    rho=0.5,k=10,eta=0.06,N=38000
  ) -> p
  plan(multisession)
  foreach (
    theta=iter(p,"row"),
    .combine=rbind,
    .options.future=list(seed=421776444)
  ) %dofuture% {
    measSIR |> pfilter(params=theta,Np=5000) -> pf
    theta$loglik <- logLik(pf)
    theta
  } -> p
  p |> arrange(Beta,mu_IR)
})-> p

p |>
  group_by(Beta,mu_IR) |>
  summarize(loglik=logmeanexp(loglik)) |>
  ungroup() |>
  mutate(loglik=ifelse(loglik>max(loglik)-25,loglik,NA)) |>
  ggplot(aes(x=Beta,y=mu_IR,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(mu[IR]))
