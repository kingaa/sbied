params <-
list(prefix = "direct/")

library(tidyverse)
library(pomp)
stopifnot(packageVersion("pomp")>="3.0")
set.seed(594709947L)

source("https://kingaa.github.io/sbied/pfilter/model.R")

measSIR %>%
  simulate(nsim=5000,format="arrays") -> x
sims <- coef(measSIR,"rho")*x$states["H",,]
matplot(time(measSIR),t(sims[1:50,]),type='l',lty=1,
  xlab="time",ylab=expression(rho*H),bty='l',col='blue')
lines(time(measSIR),obs(measSIR,"reports"),lwd=2,col='black')

ell <- dmeasure(measSIR,y=obs(measSIR),x=x$states,times=time(measSIR),log=TRUE,
  params=coef(measSIR))
dim(ell)

ell <- apply(ell,1,sum)
summary(ell)
summary(exp(ell))
