## ----prelims,include=FALSE,cache=FALSE-----------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
)

stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0.3")

set.seed(594709947L)
library(tidyverse)
library(pomp)

## ----get-data,include=FALSE----------------------------------------------
read_csv("https://kingaa.github.io/sbied/ebola/ebola_data.csv") -> dat

dat

## ----popsizes,include=FALSE----------------------------------------------
populations <- c(Guinea=10628972,Liberia=4092310,SierraLeone=6190280)

## ----plot-data,echo=FALSE------------------------------------------------
dat %>%
  ggplot(aes(x=date,y=cases,group=country,color=country))+
  geom_line()

## ----rproc,include=FALSE-------------------------------------------------
rSim <- Csnippet("
  double lambda, beta;
  double *E = &E1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Transitions
  // From class S
  double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections
  // From class E
  double transE[nstageE]; // No of transitions between classes E
  for(i = 0; i < nstageE; i++){
    transE[i] = rbinom(E[i], 1.0 - exp(-nstageE * alpha * dt));
  }
  // From class I
  double transI = rbinom(I, 1.0 - exp(-gamma * dt)); // No of transitions I->R

  // Balance the equations
  S -= transS;
  E[0] += transS - transE[0];
  for(i=1; i < nstageE; i++) {
    E[i] += transE[i-1] - transE[i];
  }
  I += transE[nstageE-1] - transI;
  R += transI;
  N_EI += transE[nstageE-1]; // No of transitions from E to I
  N_IR += transI; // No of transitions from I to R
")

rInit <- Csnippet("
  double m = N/(S_0+E_0+I_0+R_0);
  double *E = &E1;
  int j;
  S = nearbyint(m*S_0);
  for (j = 0; j < nstageE; j++) E[j] = nearbyint(m*E_0/nstageE);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);
  N_EI = 0;
  N_IR = 0;
")

## ----skel,include=FALSE--------------------------------------------------
skel <- Csnippet("
  double lambda, beta;
  const double *E = &E1;
  double *DE = &DE1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Balance the equations
  DS = - lambda * S;
  DE[0] = lambda * S - nstageE * alpha * E[0];
  for (i=1; i < nstageE; i++)
    DE[i] = nstageE * alpha * (E[i-1]-E[i]);
  DI = nstageE * alpha * E[nstageE-1] - gamma * I;
  DR = gamma * I;
  DN_EI = nstageE * alpha * E[nstageE-1];
  DN_IR = gamma * I;
")

## ----measmodel,include=FALSE---------------------------------------------
dObs <- Csnippet("
  double f;
  if (k > 0.0)
    f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
  else
    f = dpois(nearbyint(cases),rho*N_EI,1);
  lik = (give_log) ? f : exp(f);
")

rObs <- Csnippet("
  if (k > 0) {
    cases = rnbinom_mu(1.0/k,rho*N_EI);
  } else {
    cases = rpois(rho*N_EI);
  }")

## ----pomp-construction,include=FALSE-------------------------------------
ebolaModel <- function (country=c("Guinea", "SierraLeone", "Liberia"),
  timestep = 0.1, nstageE = 3) {

  ctry <- match.arg(country)
  pop <- unname(populations[ctry])
  nstageE <- as.integer(nstageE)

  globs <- paste0("static int nstageE = ",nstageE,";")

  dat %>% filter(country==ctry) %>% select(-country) -> dat

  ## Create the pomp object
  dat %>%
    select(week,cases) %>%
    pomp(
      times="week",
      t0=min(dat$week)-1,
      globals=globs,
      accumvars=c("N_EI","N_IR"),
      statenames=c("S",sprintf("E%1d",seq_len(nstageE)),
        "I","R","N_EI","N_IR"),
      paramnames=c("N","R0","alpha","gamma","rho","k",
        "S_0","E_0","I_0","R_0"),
      dmeasure=dObs, rmeasure=rObs,
      rprocess=discrete_time(step.fun=rSim, delta.t=timestep),
      skeleton=vectorfield(skel),
      partrans=parameter_trans(
        log=c("R0","k"),logit="rho",
        barycentric=c("S_0","E_0","I_0","R_0")),
      rinit=rInit
    ) -> po
}

ebolaModel("Guinea") -> gin
ebolaModel("SierraLeone") -> sle
ebolaModel("Liberia") -> lbr

## ----load-profile,echo=FALSE---------------------------------------------
read_csv("https://kingaa.github.io/sbied/ebola/ebola_profiles.csv") -> profs

## ----profiles-plots,results='hide',echo=FALSE----------------------------
library(tidyverse)
theme_set(theme_bw())

profs %>%
  gather(variable,value,-profile,-country,-loglik) %>%
  filter(variable==profile) %>%
  group_by(country) %>%
  mutate(dll=loglik-max(loglik)) %>%
  group_by(country,profile,value) %>%
  filter(loglik==max(loglik)) %>%
  ungroup() %>%
  ggplot(aes(x=value,y=dll))+
  geom_point(color="red")+
  geom_hline(yintercept=-0.5*qchisq(p=0.99,df=1))+
  facet_grid(country~profile,scales="free")+
  labs(y=expression(l))

## ----diagnostics1a,echo=FALSE---------------------------------------------
profs %>%
  filter(country=="Guinea") %>%
  filter(loglik==max(loglik)) %>%
  select(-loglik,-loglik.se,-country,-profile) -> coef(gin)

## ----diagnostics1b,echo=FALSE---------------------------------------------
gin %>%
  simulate(nsim=20,format="data.frame",include.data=TRUE) %>%
  mutate(
    date=min(dat$date)+7*(week-1),
    is.data=ifelse(.id=="data","yes","no")
  ) %>%
  ggplot(aes(x=date,y=cases,group=.id,color=is.data,alpha=is.data))+
  geom_line()+
  guides(color=FALSE,alpha=FALSE)+
  scale_color_manual(values=c(no=gray(0.6),yes="red"))+
  scale_alpha_manual(values=c(no=0.5,yes=1))

## ----diagnostics-growth-rate---------------------------------------------
growth.rate <- function (y) {
  cases <- y["cases",]
  fit <- lm(log1p(cases)~seq_along(cases))
  unname(coef(fit)[2])
}

gin %>%
  probe(probes=list(r=growth.rate),nsim=500) %>%
  plot()

## ----diagnostics-growth-rate-and-sd--------------------------------------
growth.rate.plus <- function (y) {
  cases <- y["cases",]
  fit <- lm(log1p(cases)~seq_along(cases))
  c(r=unname(coef(fit)[2]),sd=sd(residuals(fit)))
}

gin %>%
  probe(probes=list(growth.rate.plus),nsim=500) %>%
  plot()

## ----diagnostics2,fig.height=6-------------------------------------------
log1p.detrend <- function (y) {
  cases <- y["cases",]
  fit <- lm(log1p(cases)~seq_along(cases))
  y["cases",] <- as.numeric(residuals(fit))
  y
}

gin %>%
  probe(nsim=500,
    probes=list(
      growth.rate.plus,
      probe.quantile(var="cases",prob=c(0.25,0.75)),
      probe.acf(var="cases",lags=c(1,2),type="correlation",
        transform=log1p.detrend))) %>%
  plot()

## ----forecasts1a----------------------------------------------------------
library(pomp)
library(tidyverse)

set.seed(988077383L)

## forecast horizon
horizon <- 13

## ----forecasts1b----------------------------------------------------------

## Weighted quantile function
wquant <- function (x, weights, probs = c(0.025,0.5,0.975))
{
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  w <- cumsum(weights)/sum(weights)
  rval <- approx(w,x,probs,rule=1)
  rval$y
}

## ----forecasts1c----------------------------------------------------------
profs %>%
  filter(country=="SierraLeone") %>%
  select(-country,-profile,-loglik.se) %>%
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.99)) %>%
  gather(parameter,value) %>%
  group_by(parameter) %>%
  summarize(min=min(value),max=max(value)) %>%
  ungroup() %>%
  filter(parameter!="loglik") %>%
  column_to_rownames("parameter") %>%
  as.matrix() -> ranges

## ----forecasts1d----------------------------------------------------------
sobol_design(
  lower=ranges[,"min"],
  upper=ranges[,"max"],
  nseq=20
) -> params
plot(params)

## ----forecasts2a----------------------------------------------------------
bake(file="forecasts.rds",{
library(foreach)
library(doParallel)
library(iterators)
library(doRNG)

registerDoParallel()
registerDoRNG(887851050L)

## ----forecasts2b----------------------------------------------------------
foreach(p=iter(params,by="row"),
  .inorder=FALSE,
  .combine=bind_rows
) %dopar% {

  library(pomp)

## ----forecasts2c----------------------------------------------------------
  M1 <- ebolaModel("SierraLeone")

  M1 %>% pfilter(params=p,Np=2000,save.states=TRUE) -> pf

## ----forecasts2d----------------------------------------------------------
  pf %>%
    saved.states() %>% ## latent state for each particle
    tail(1) %>%        ## last timepoint only
    melt() %>%         ## reshape and rename the state variables
    spread(variable,value) %>%
    group_by(rep) %>%
    summarize(S_0=S, E_0=E1+E2+E3, I_0=I, R_0=R) %>%
    gather(variable,value,-rep) %>%
    spread(rep,value) %>%
    column_to_rownames("variable") %>%
    as.matrix() -> x

## ----forecasts2e1----------------------------------------------------------
  pp <- parmat(unlist(p),ncol(x))

## ----forecasts2e2----------------------------------------------------------
  M1 %>%
    simulate(params=pp,format="data.frame") %>%
    select(.id,week,cases) %>%
    mutate(
      period="calibration",
      loglik=logLik(pf)
    ) -> calib

## ----forecasts2f----------------------------------------------------------
  M2 <- M1
  time(M2) <- max(time(M1))+seq_len(horizon)
  timezero(M2) <- max(time(M1))

## ----forecasts2g----------------------------------------------------------
  pp[rownames(x),] <- x

  M2 %>%
    simulate(params=pp,format="data.frame") %>%
    select(.id,week,cases) %>%
    mutate(
      period="projection",
      loglik=logLik(pf)
    ) -> proj

## ----forecasts2h----------------------------------------------------------
  bind_rows(calib,proj) -> sims

## ----forecasts2i----------------------------------------------------------
}}) -> sims

## ----forecasts2j----------------------------------------------------------
sims %>%
  mutate(weight=exp(loglik-mean(loglik))) %>%
  arrange(week,.id) -> sims

## ----forecasts2k----------------------------------------------------------
sims %>%
  filter(week==max(week)) %>%
  summarize(ess=sum(weight)^2/sum(weight^2))

## ----forecasts2l----------------------------------------------------------
sims %>%
  group_by(week,period) %>%
  summarize(
    p=c(0.025,0.5,0.975),
    q=wquant(cases,weights=weight,probs=p),
    label=c("lower","median","upper")
  ) %>%
  select(-p) %>%
  spread(label,q) %>%
  ungroup() %>%
  mutate(date=min(dat$date)+7*(week-1)) -> simq

## ----forecast-plots,echo=FALSE-------------------------------------------
simq %>%
  ggplot(aes(x=date))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=period),alpha=0.3,color=NA)+
  geom_line(aes(y=median,color=period))+
  geom_point(data=filter(dat,country=="SierraLeone"),
    mapping=aes(x=date,y=cases),color="black")+
  labs(y="cases")
