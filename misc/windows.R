params <-
list(prefix = "windows")

set.seed(1350254336)

options(pomp_cdir="./tmp")

library(tidyverse)
library(pomp)

library(doParallel)
d <- detectCores()
registerDoParallel(d-1)

source("https://kingaa.github.io/sbied/pfilter/model.R")

foreach (i=1:4,.combine=c) %dopar% {
  library(pomp)
  measSIR %>%
    mif2(
      Np=1000, Nmif=5,
      cooling.fraction.50=0.5,
      rw.sd=rw.sd(Beta=0.2, rho=0.2, eta=ivp(0.2)),
### Compilation is triggered here, by the call to `parameter_trans`.
      partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
      paramnames=c("Beta","rho","eta")
    )
} -> mifs_local

measSIR %>%
  pomp(
### Compilation is triggered here, outside the parallel block.
    partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
    paramnames=c("Beta","rho","eta")
  ) -> measSIR2

foreach (i=1:4,.combine=c) %dopar% {
### No compilation is triggered inside the parallel code block.
  library(pomp)
  measSIR2 %>%
    mif2(
      Np=1000, Nmif=5,
      cooling.fraction.50=0.5,
      rw.sd=rw.sd(Beta=0.2, rho=0.2, eta=ivp(0.2))
    )
} -> mifs_local
