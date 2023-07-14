params <-
list(prefix = "windows")

set.seed(1350254336)

options(pomp_cdir="./tmp")

library(tidyverse)
library(pomp)
library(doFuture)

source("https://kingaa.github.io/sbied/pfilter/model.R")

plan(multisession)
foreach (
  i=1:4,
  .combine=c,
  .options.future=list(seed=TRUE)
) %dofuture% {
  measSIR |>
    mif2(
      Np=1000, Nmif=5,
      cooling.fraction.50=0.5,
      rw.sd=rw_sd(Beta=0.2, rho=0.2, eta=ivp(0.2)),
### Compilation is triggered here, by the call to `parameter_trans`.
      partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
      paramnames=c("Beta","rho","eta")
    )
} -> mifs_local

measSIR |>
  pomp(
### Compilation is triggered here, outside the parallel block.
    partrans=parameter_trans(log="Beta",logit=c("rho","eta")),
    paramnames=c("Beta","rho","eta")
  ) -> measSIR2

plan(multisession)
foreach (
  i=1:4,
  .combine=c,
  .options.future=list(seed=TRUE)
) %dofuture% {
### No compilation is triggered inside the parallel code block.
  measSIR2 |>
    mif2(
      Np=1000, Nmif=5,
      cooling.fraction.50=0.5,
      rw.sd=rw_sd(Beta=0.2, rho=0.2, eta=ivp(0.2))
    )
} -> mifs_local
