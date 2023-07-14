library(doFuture)
plan(multisession)
set.seed(2488820)

library(tidyverse)
read_csv(paste0("https://kingaa.github.io/sbied/stochsim/",
  "Measles_Consett_1948.csv")) |>
  select(week,reports=cases) -> meas
meas |> as.data.frame() |> head()

library(tidyverse)
meas |>
  ggplot(aes(x=week,y=reports))+
  geom_line()+
  geom_point()

sir_step <- function (S, I, R, N, Beta, mu_IR, delta.t, ...)
{
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  c(S = S, I = I, R = R)
}

sir_rinit <- function (N, eta, ...) {
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)))
}

library(pomp)
meas |>
  pomp(times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_rinit
  ) -> measSIR

sir_step <- function (S, I, R, N, Beta, mu_IR, delta.t,
  H, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR
  c(S = S, I = I, R = R, H = H)
}

sir_rinit <- function (N, eta, ...) {
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)), H = 0)
}

measSIR |>
  pomp(
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_rinit, accumvars="H"
  ) -> measSIR

sir_dmeas <- function (reports, H, rho, k, log, ...) {
  dnbinom(x=reports, size=k, mu=rho*H, log=log)
}

sir_rmeas <- function (H, rho, k, ...) {
  c(reports=rnbinom(n=1, size=k, mu=rho*H))
}

measSIR |>
  pomp(
    rmeasure=sir_rmeas,
    dmeasure=sir_dmeas
  ) -> measSIR



sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_rinit <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

sir_dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,rho*H,give_log);
")

sir_rmeas <- Csnippet("
  reports = rnbinom_mu(k,rho*H);
")

measSIR |>
  pomp(rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_rinit,
    rmeasure=sir_rmeas,
    dmeasure=sir_dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","mu_IR","N","eta","rho","k")
  ) -> measSIR

measSIR |>
  simulate(
    params=c(Beta=7.5,mu_IR=0.5,rho=0.5,k=10,
      eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE
  ) -> sims

sims |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")
