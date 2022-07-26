params <-
list(prefix = "pfnm")

library(tidyverse)
library(pomp)
set.seed(594709947L)

source("https://kingaa.github.io/sbied/pfilter/model.R")

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


stew(file="like_optim2.rda",{
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
})

mle %>% simulate(nsim=10,format="data.frame",include.data=TRUE) -> sims

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  guides(color="none")+
  geom_line()
