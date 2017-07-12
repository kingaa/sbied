#' ---
#' title: "Overdispersion on the boarding school flu data"
#' author: "Edward Ionides"
#' output:
#'   html_document:
#'     toc: yes
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
#' ---
#' 
#' \newcommand\prob[1]{\mathbb{P}\left[{#1}\right]}
#' \newcommand\expect[1]{\mathbb{E}\left[{#1}\right]}
#' \newcommand\var[1]{\mathrm{Var}\left[{#1}\right]}
#' \newcommand\dist[2]{\mathrm{#1}\left(#2\right)}
#' \newcommand\dlta[1]{{\Delta}{#1}}
#' \newcommand\lik{\mathcal{L}}
#' \newcommand\loglik{\ell}
#' 
#' [Licensed under the Creative Commons Attribution-NonCommercial license](http://creativecommons.org/licenses/by-nc/4.0/).
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' Produced with **R** version `r getRversion()` and **pomp** version `r packageVersion("pomp")`.
#' 
## ----prelims,include=FALSE,purl=TRUE,cache=FALSE-------------------------
library(pomp)
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="1.12")
set.seed(557976883)

#' 
#' -----------
#' 
#' ----------
#' 
#' ## Introduction
#' 
#' - This tutorial investigates the consequences of modeling dynamic overdispersion on the boarding school flu dataset.
#' 
#' - It extends the example applying `mif` to the flu data in the [iterated filtering tutorial](../mif/mif.html).
#' 
#' <br>
#' 
#' -----
#' 
#' ----
#' 
#' ## Objectives
#' 
#' 1. Provide a relatively simple example to better understand overdispersed Markov chain models that are essential for large populations such as the [city-level measles case study](../measles/measles.html).
#' 
#' 2. Investigate whether fitting models with dynamic overdispersion (i.e., environmental stochasticity, also known as extra-demographic stochasticity) is helpful even for this relatively small population, for which demographic stochasticity is relatively large.
#' 
#' 3. Consider the respective roles of overdispersion in the measurement model and the dynamic model.
#' 
#' 4. Think about scientific interpretations of overdispersion if it is statistically evident.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ## Adding overdispersion to the basic SIR model
#' 
#' 
#' Recall the data on an influenza outbreak in a British boarding school [@anonymous78].
#' Reports consist of the number of children confined to bed for each of the 14 days of the outbreak.
#' The total number of children at the school was 763, and a total of 512 children spent time away from class.
#' Only one adult developed influenza-like illness, so adults are omitted from the data and model.
#' First, we read in the data:
#' 
## ----load_bbs------------------------------------------------------------
bsflu_data <- read.table("bsflu_data.txt")

#' 
#' Our model is a variation on a basic SIR Markov chain, with state $X(t)=(S(t),I(t),R_1(t),R_2(t),R_3(t))$ giving the numbers of individuals in the susceptible and infectious categories, and three stages of recovery.
#' The recovery stages, $R_1$, $R_2$ and $R_3$, are all modeled to be non-contagious.
#' $R_1$ consists of individuals who are bed-confined if they showed symptoms;
#' $R_2$ consists of individuals who are convalescent if they showed symptoms;
#' $R_3$ consists of recovered individuals who have returned to school-work if they were symtomatic.
#' The observation on day $n$ of the observed epidemic (with $t_1$ being 22 January) consists of the numbers of children who are bed-confined and convalescent.
#' Ten individuals received antibiotics for secondary infections, and they had longer bed-confinement and convalescence times.
#' Partly for this reason, and because our primary interest is in parameters related to transmission, we'll narrow our focus to the bed-confinement numbers, $B_n$, modeling these as $B_n\sim\dist{Poisson}{\rho R_1(t_n)}$, where $\rho$ is a reporting rate corresponding to the chance an infected boy is symptomatic.
#' 
#' 
#' The index case for the epidemic was purportedly a boy recently returned from Hong Kong, who was reported to have a transient febrile illness from 15 to 18 January.
#' It would therefore be reasonable to initialize the epidemic at $t_0=-6$ with $I(t_0)=1$.
#' This is a little tricky to reconcile with the rest of the data;
#' for the moment, we avoid this issue by instead initializing with $I(t_0)=1$ at $t_0=0$.
#' All other individuals are modeled to be initially susceptible.
#' 
#' Our previous Markov transmission model is that each individual in $S$ transitions to $I$ at rate $\beta\,I(t)/N$.
#' We are going to extend this to include stochastic variation by incorporating multiplicative noise,
#' $$\mu_{SI}=\beta\frac{I(t)}{N}d\Gamma/dt.$$
#' 
#' * Here, $\Gamma(t)$ is a gamma process with $\E[\Gamma(t)]=t$ and $\var[\Gamma(t)]=\sigma^2 t$.
#' 
#' Thus, $\sigma^2$ is the __infinitesimal variance parameter__ of the noise process, which we will call the __extrademographic process noise__ parameter.
#' 
#' * Why do we use multiplicative noise, rather than additive noise such as 
#' $$\mu_{SI}^\prime=\beta\frac{I(t)}{N} + d\Gamma/dt?$$
#' 
#' * The gamma process has the property of being non-negative, which is the main reason it may be preferable to Gaussian noise.
#' 
#' * The gamma process is a pure jump process, constant between jumps. There are an infinite number of jumps in any finite time interval, but almost all of them are negligibly small. Thus, the derivative of the gamma process doesn't exist in the usual sense, but one can still give formal meaning to $d\Gamma/dt$. All this is similar to Gaussian noise, which can also be considered as the formal derivative of a non-differentiable process (Brownian motion).
#' 
#' Over-dispersion in the measurement model may be appropriate for similar reasons. Measurement over-dispersion might be considered together with, or in place of, dynamic over-dispersion. 
#' 
#' * Adding gamma noise to a Poisson measurement model leads to a negative binomial measurement model. The overdispersion parameter is $\psi$, with the variance for mean $\mu$ being $\mu+\mu^2/\psi$ and the Poisson model being recovered in the limit as $\psi\to\infty$.
#' 
#' The remainder of the model is unchanged, with each individual in $I$ transitioning at rate $\mu_I$ to $R_1$, subsequently moving from $R_1$ to $R_2$ at  rate $\mu_{R_1}$, and finally from $R_2$ to $R_3$ at rate $\mu_{R_2}$.
#' Therefore, $1/\mu_I$ is the mean infectious time prior to bed-confinement; $1/\mu_{R_1}$ is the mean duration of bed-confinement for symptomatic cases;
#' $1/\mu_{R_2}$ is the mean duration of convalescence for symptomatic cases.
#' All rates have units $\mathrm{day}^{-1}$. 
#' 
#' We do not need a representation of $R_3$ since this variable has consequences neither for the dynamics of the state process nor for the data.
#' Since we are confining ourselves for the present to fitting only the $B_n$ data, we need not track $R_2$.
#' We enumerate the state variables ($S$, $I$, $R_1$) and the parameters ($\beta$, $\mu_I$, $\rho$, $\mu_{R_1}$, $\sigma$, $\psi$) as follows:
#' 
## ----bsflu_names---------------------------------------------------------
statenames <- c("S","I","R1")
paramnames <- c("Beta","mu_I","mu_R1","rho","sigma","psi")

#' 
#' In the codes below, we'll refer to the data variables by their names ($B$, $C$), as given in the `bsflu_data` data-frame:
#' 
#' Now, we write the model code:
#' 
## ----csnippets_bsflu-----------------------------------------------------
### for comparison, here is the redundant equidispersed model
dmeas <- Csnippet("
  lik = dpois(B,rho*R1+1e-6,give_log);
")
rmeas <- Csnippet("
  B = rpois(rho*R1+1e-6);
")
########

dmeas <- Csnippet("
  lik = dnbinom(B, psi, psi/(psi+rho*R1+1e-6), give_log);
")

rmeas <- Csnippet("
  B = rnbinom(psi, psi/(psi+rho*R1+1e-6));
")

rproc <- Csnippet("
  double dGamma;
  double N = 763;
  dGamma = rgammawn(sigma,dt);
  double t1 = rbinom(S,1-exp(-Beta*I/N*dGamma));
  double t2 = rbinom(I,1-exp(-mu_I*dt));
  double t3 = rbinom(R1,1-exp(-mu_R1*dt));
  S  -= t1;
  I  += t1 - t2;
  R1 += t2 - t3;
")

init <- Csnippet("
 S = 762;
 I = 1;
 R1 = 0;
")

fromEst <- Csnippet("
 TBeta = exp(Beta);
 Tmu_I = exp(mu_I);
 Trho = expit(rho);
 Tsigma = exp(sigma);
 Tpsi = exp(psi);
")

toEst <- Csnippet("
 TBeta = log(Beta);
 Tmu_I = log(mu_I);
 Trho = logit(rho);
 Tsigma = log(sigma);
 Tpsi = log(psi);
")

#' 
#' Now we build the `pomp` object:
#' 
## ----pomp_bsflu----------------------------------------------------------
library(pomp)

pomp(
  data=subset(bsflu_data,select=-C),
  times="day",t0=0,
  rmeasure=rmeas,dmeasure=dmeas,
  rprocess=euler.sim(rproc,delta.t=1/12),
  initializer=init,
  fromEstimationScale=fromEst,toEstimationScale=toEst,
  statenames=statenames,
  paramnames=paramnames
) -> bsflu

#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Testing the codes.
#' 
#' To develop and debug code, it is useful to have testing codes that run quickly and fail if the codes are not working correctly.
#' As such a test, here we run some simulations and a particle filter.
#' We'll use the following parameters, derived from our earlier explorations:
## ----start_params--------------------------------------------------------
params <- c(Beta=2,mu_I=1,rho=0.9,mu_R1=1/3,mu_R2=1/2,sigma=0.2,psi=20)

#' 
#' Now to run and plot some simulations:
## ----init_sim------------------------------------------------------------
y <- simulate(bsflu,params=params,nsim=10,as.data.frame=TRUE)

#' 
#' Before engaging in iterated filtering, it is a good idea to check that the basic particle filter is working since iterated filtering builds on this technique.
#' The simulations above check the `rprocess` and `rmeasure` codes;
#' the particle filter depends on the `rprocess` and `dmeasure` codes and so is a check of the latter.
#' 
#' We need to find, by trial and error, a suitable number of particles making a compromise between (i) accuracy, which depends also on the model and data; (ii) computational resources; (iii) the timescale on which we want an answer. 
## ----np------------------------------------------------------------------
NP <- 2000

#' It is helpful for debugging to have a flag that switches the algorithmic parameters to makes the code run quickly.
#' For code development, it is also helpful to have a flag that gives algorithmic settings which run in a convenient time, say 30 minutes.
## ----debug-np, include=TRUE----------------------------------------------
SHORT_RUN <- TRUE
if(SHORT_RUN) NP <- 1000
#DEBUG <- TRUE
DEBUG <- FALSE
if(DEBUG) NP <- 50

#' 
#' We will not report the debugging values of subsequent parameters.
#' Now we compute the likelihood at our parameter guess.
## ----init_pfilter--------------------------------------------------------
pf <- pfilter(bsflu,params=params,Np=NP)

#' 
#' The above plot shows the data (`B`), along with the *effective sample size* of the particle filter (`ess`) and the log likelihood of each observation conditional on the preceding ones (`cond.logLik`).
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Setting up the estimation problem.
#' 
#' Let's treat $\mu_{R_1}$ and  $\mu_{R_2}$ as known, and fix these parameters at the empirical means of the bed-confinement and convalescence times for symptomatic cases, respectively:
#' 
## ----fixed_params--------------------------------------------------------
(fixed_params <- with(bsflu_data,c(mu_R1=1/(sum(B)/512),mu_R2=1/(sum(C)/512))))

#' 
#' We will estimate $\beta$, $\mu_I$, $\rho$ and $\sigma$.
#' 
#' It will be helpful to parallelize most of the computations.
#' Most machines nowadays have multiple cores and using this computational capacity is as simple as:
#' 
#' i. letting **R** know you plan to use multiple processors;
#' i. using the parallel for loop provided by the **foreach** package; and
#' i. paying proper attention to the use of parallel random number generators.
#' 
#' For example:
#' 
## ----parallel-setup,cache=FALSE------------------------------------------
library(foreach)
library(doParallel)
registerDoParallel()

#' 
#' The first two lines above load the **foreach** and **doParallel** packages, the latter being a "backend" for the **foreach** package.
#' The next line tells **foreach** that we will use the **doParallel** backend.
#' By default, **R** will guess how many cores are available and will run about half this number of concurrent **R** processes.
#' 
#' ### Running a particle filter.
#' 
#' We proceed to carry out replicated particle filters at an initial guess of $\beta=2$, $\mu_I=1$, $\rho=0.9$ and $\sigma=0.2$.
#' 
#' 
## ----pf------------------------------------------------------------------
library(doRNG)
registerDoRNG(625904618)
bake(file="pf.rds",{
  foreach(i=1:10,.packages='pomp',
    .export=c("bsflu","fixed_params")
  ) %dopar% {
    pfilter(bsflu,params=c(Beta=2,mu_I=1,rho=0.9,sigma=0.2,psi=20,fixed_params),Np=NP)
  }
}) -> pf
(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))

#' 
#' In `r round(attr(pf,"system.time")["elapsed"],2)` seconds, using `r min(getDoParWorkers(),length(pf))` cores, we obtain an unbiased likelihood estimate of `r round(L_pf[1],1)` with a Monte Carlo standard error of `r signif(L_pf[2],2)`.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### A local search of the likelihood surface
#' 
#' Let's carry out a local search using `mif2` around this point in parameter space. 
#' To do so, we need to choose the `rw.sd` and `cooling.fraction.50` algorithmic parameters.
#' Since $\beta$ and $\mu_I$ will be estimated on the log scale, and we expect that multiplicative perturbations of these parameters will have roughly similar effects on the likelihood, we'll use a perturbation size of $0.02$, which we imagine will have a small but non-negligible effect.
#' For simplicity, we'll use the same perturbation size on $\rho$.
#' We fix `cooling.fraction.50=0.5`, so that after 50 `mif2` iterations, the perturbations are reduced to half their original magnitudes.
#' 
## ----nmif----------------------------------------------------------------
NMIF <- 100
NP_MIF <- NP/2

## ----debug-nmif, include=FALSE-------------------------------------------
if(SHORT_RUN){
  NP_MIF <- NP/2
  NMIF <- 50
}
if(DEBUG) {
  NP_MIF <- 50
  NMIF <- 5
}

#' 
## ----box_search_local----------------------------------------------------
registerDoRNG(482947940)
bake(file="box_search_local.rds",{
foreach(i=1:20,
    .packages='pomp',
    .combine=c, 
    .export=c("bsflu","fixed_params")
  ) %dopar% {
    mif2(
      bsflu,
      start=c(Beta=2,mu_I=1,rho=0.9,sigma=0.2,psi=20,fixed_params),
      Np=NP_MIF,
      Nmif=NMIF,
      cooling.type="geometric",
      cooling.fraction.50=0.5,
      transform=TRUE,
      rw.sd=rw.sd(Beta=0.02,mu_I=0.02,rho=0.02,sigma=0.02,psi=0.02)
    )
  } 
}) -> mifs_local

#' 
#' We obtain some diagnostic plots with the `plot` command applied to `mifs_local`.
#' Here is a way to get a prettier version:
#' 
#' 
#' No filtering failures (`nfail`) are generated at any point, which is comforting.
#' In general, we expect to see filtering failures whenever our initial guess (`start`) is incompatible with one or more of the observations.
#' Filtering failures at the MLE are an indication that the model, at its best, is incompatible with one or more of the data.
#' 
#' We see that the likelihood generally increases as the iterations proceed, though there is considerable variability due to the stochastic nature of this Monte Carlo algorithm.
#' Although the filtering carried out by `mif2` in the final filtering iteration generates an approximation to the likelihood at the resulting point estimate, this is not usually good enough for reliable inference.
#' Partly, this is because parameter perturbations are applied in the last filtering iteration, so that the likelhood shown here is not identical to that of the model of interest.
#' Partly, this is because `mif2` is usually carried out with fewer particles than are needed for a good likelihood evaluation:
#' the errors in `mif2` average out over many iterations of the filtering.
#' Therefore, we evaluate the likelihood, together with a standard error, using replicated particle filters at each point estimate:
#' 
## ----lik_local-----------------------------------------------------------
registerDoRNG(900242057)
bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.packages='pomp',.combine=rbind) %dopar% 
  {
    evals <- replicate(10, logLik(pfilter(mf,Np=NP)))
    ll <- logmeanexp(evals,se=TRUE)
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }
}) -> results_local

## ------------------------------------------------------------------------
results_local <- as.data.frame(results_local)

#' 
#' This investigation took  `r round(attr(mifs_local,"system.time")["elapsed"],0)` sec for the maximization and `r round(t_local["elapsed"],0)` sec for the likelihood evaluation.
#' These repeated stochastic maximizations can also show us the geometry of the likelihood surface in a neighborhood of this point estimate:
#' 
#' 
#' Although this plot some hints of ridges in the likelihood surface (cf. the $\beta$-$\mu_I$ panel), the sampling is still too sparse to give a clear picture.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### A global search of the likelihood surface using randomized starting values
#' 
#' When carrying out parameter estimation for dynamic systems, we need to specify beginning values for both the dynamic system (in the state space) and the parameters (in the parameter space).
#' To avoid confusion, we use the term "initial values" to refer to the state of the system at $t_0$ and "starting values" to refer to the point in parameter space at which a search is initialized.
#' 
#' Practical parameter estimation involves trying many starting values for the parameters.
#' One way to approach this is to choose a large box in parameter space that contains all remotely sensible parameter vectors.
#' If an estimation method gives stable conclusions with starting values drawn randomly from this box, this gives some confidence that an adequate global search has been carried out. 
#' 
#' For our flu model, a box containing reasonable parameter values might be
#' 
## ----box_global----------------------------------------------------------
params_box <- rbind(
  Beta=c(1,5),
  mu_I=c(0.5,3),
  rho = c(0.5,1),
  sigma = c(0.01,0.2),
  psi = c(1,20)
)

#' 
#' We are now ready to carry out likelihood maximizations from
## ----nglobal-------------------------------------------------------------
NGLOBAL <- 300

#' diverse starting points.
#' 
## ----debug-nglobal, include=FALSE----------------------------------------
if(SHORT_RUN){
  NGLOBAL <- 100
}
if(DEBUG) {
  NGLOBAL <- 10
}

#' 
## ----box_search_global---------------------------------------------------
registerDoRNG(1270401374)
guesses <- as.data.frame(apply(params_box,1,function(x)runif(NGLOBAL,x[1],x[2])))
mf1 <- mifs_local[[1]]
bake(file="box_search_global.rds",{
  foreach(guess=iter(guesses,"row"), 
    .packages='pomp', 
    .combine=rbind,
    .options.multicore=list(set.seed=TRUE),
    .export=c("mf1","fixed_params")
  ) %dopar% 
  {
    mf <- mif2(mf1,start=c(unlist(guess),fixed_params))
    mf <- mif2(mf,Nmif=NMIF)
    ll <- replicate(10,logLik(pfilter(mf,Np=NP_MIF)))
    ll <- logmeanexp(ll,se=TRUE)
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }
}) -> results_global

## ------------------------------------------------------------------------
results_global <- as.data.frame(results_global)

#' The above codes run one search from each of `r nrow(guesses)` starting values.
#' Each search consissts of an initial run of `r nrow(conv.rec(mf1))` IF2 iterations, followed by another 100 iterations.
#' These codes exhibit a general **pomp** behavior:
#' re-running a command on an object (i.e., `mif2` on `mf1`) created by the same command preserves the algorithmic arguments.
#' In particular, running `mif2` on the result of a `mif2` computation re-runs IF2 from the endpoint of the first run.
#' In the second computation, by default, all algorithmic parameters are preserved;
#' here we overrode the default choice of `Nmif`.
#' 
#' Following the `mif2` computations, the particle filter is used to evaluate the likelihood, as before.
#' In contract to the local-search codes above, here we return only the endpoint of the search, together with the likelihood estimate and its standard error in a named vector.
#' The best result of this search had a likelihood of `r round(max(results_global$loglik),1)` with a standard error of `r round(results_global$loglik.se[which.max(results_global$loglik)],2)`.
#' This took `r round(t_global["elapsed"]/60,1)` minutes altogether using `r n_global` processors.
#' 
#' Again, we attempt to visualize the global geometry of the likelihood surface using a scatterplot matrix.
#' In particular, here we plot both the starting values (grey) and the IF2 estimates (red).
#' 
#' 
#' We see that optimization attempts from diverse remote starting points converge on a particular region in parameter space.
#' Moreover, the estimates have comparable likelihoods, despite their considerable variability.
#' This gives us some confidence in our maximization procedure. 
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Do we need both dynamic noise and measurement noise?
#' 
#' Let's try with only $\sigma$.
## ----box_search_global_sigma---------------------------------------------
params_box_sigma <- params_box[c("Beta","mu_I","rho","sigma"),]
fixed_params_sigma <- c(fixed_params,psi=1e6)
registerDoRNG(1270401374)
guesses <- as.data.frame(apply(params_box_sigma,1,function(x)runif(NGLOBAL,x[1],x[2])))
mf1 <- mifs_local[[1]]
bake(file="box_search_global_sigma.rds",{
  foreach(guess=iter(guesses,"row"), 
    .packages='pomp', 
    .combine=rbind,
    .options.multicore=list(set.seed=TRUE),
    .export=c("mf1","fixed_params_sigma")
  ) %dopar% 
  {
    mf <- mif2(mf1,
        start=c(unlist(guess),fixed_params_sigma),
	Nmif=NMIF,
	rw.sd=rw.sd(Beta=0.02,mu_I=0.02,rho=0.02,sigma=0.02)
    )
    ll <- replicate(10,logLik(pfilter(mf,Np=NP_MIF)))
    ll <- logmeanexp(ll,se=TRUE)
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }
}) -> results_global_sigma

#' 
## ---- include=FALSE------------------------------------------------------
results_global_sigma <- as.data.frame(results_global_sigma)

#' 
#' The best result of this search had a likelihood of `r round(max(results_global_sigma$loglik),1)` with a standard error of `r round(results_global_sigma$loglik.se[which.max(results_global_sigma$loglik)],2)`. We see that dynamic overdispersion, by itself, gives negligible improvement in model fit.
#' 
#' Now let's try with only $\psi$.
## ----box_search_global_psi-----------------------------------------------
params_box_psi <- params_box[c("Beta","mu_I","rho","psi"),]
fixed_params_psi <- c(fixed_params,sigma=1e-6)
registerDoRNG(1270401374)
guesses <- as.data.frame(apply(params_box_psi,1,function(x)runif(NGLOBAL,x[1],x[2])))
mf1 <- mifs_local[[1]]
bake(file="box_search_global_psi.rds",{
  foreach(guess=iter(guesses,"row"), 
    .packages='pomp', 
    .combine=rbind,
    .options.multicore=list(set.seed=TRUE),
    .export=c("mf1","fixed_params_psi")
  ) %dopar% 
  {
    mf <- mif2(mf1,
        start=c(unlist(guess),fixed_params_psi),
	Nmif=NMIF,
	rw.sd=rw.sd(Beta=0.02,mu_I=0.02,rho=0.02,psi=0.02)
    )
    ll <- replicate(10,logLik(pfilter(mf,Np=NP_MIF)))
    ll <- logmeanexp(ll,se=TRUE)
    c(coef(mf),loglik=ll[1],loglik=ll[2])
  }
}) -> results_global_psi

#' 
## ---- include=FALSE------------------------------------------------------
results_global_psi <- as.data.frame(results_global_psi)

#' 
#' The best result of this search had a likelihood of `r round(max(results_global_psi$loglik),1)` with a standard error of `r round(results_global_psi$loglik.se[which.max(results_global_psi$loglik)],2)`.
#' So, measurement overdispersion substantially improves the fit, and also greatly reduces the filtering error.
#' 
#' To understand what is going on, let's look at the fitted MLE models under these different constraints
## ----MLEs----------------------------------------------------------------
mle_psi <- unlist(results_global_psi[which.max(results_global_psi$loglik),paramnames])
mle_sigma <- unlist(results_global_sigma[which.max(results_global_sigma$loglik),paramnames])
cbind(mle_psi,mle_sigma)

#' Call $H_\sigma$ the process overdispersion model and $H_\psi$ the observation overdispersion model.
#' We see that the main difference is a lower reporting rate estimate for $H_\psi$.
#' Let's see which data points are primarily responsible for the improvement in the likelihood
## ----fitted-models-------------------------------------------------------
pf_psi <- pfilter(bsflu,params=mle_psi,Np=5*NP,pred.mean=TRUE)
pf_sigma <- pfilter(bsflu,params=mle_sigma,Np=5*NP,pred.mean=TRUE)
plot(cond.logLik(pf_psi)-cond.logLik(pf_sigma))

#' 
#' We see that the big gain for $H_\psi$ is at the end of the time series.
#' Around the main peak of incidence, $H_\sigma$ actually fits better.
#' To better understand what is going on, let's look at the one-step prediction means. The differences between these and the observations are the residuals.
#' 
## ----predictions---------------------------------------------------------
plot(obs(bsflu)["B",])
lines(pred.mean(pf_psi)["R1",]*mle_psi["rho"],lty="dashed",col="blue")
lines(pred.mean(pf_sigma)["R1",]*mle_sigma["rho"],lty="dotted",col="red")


#' 
#' We see that the MLE found for $H_\sigma$ cannot explain how the cases drop off so fast, particularly on day 12.
#' The fit for $H_\psi$ may appear superficially worse, but is objectively better by the standards of likelihood.
#' We should not be entirely satisfied with either model $H_\psi$ or $H_\sigma$; a better model might be able to capture the peak as well as the early and late epidemic.
#' However, the late epidemic is perhaps of less interest, since the initial transmission and peak incidence of the epidemic may be of primary interest for disease control.
#' 
#' 
#' --------------------------
#' 
#' 
#' ## [Back to course homepage](../index.html)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/mif/mif.R)
#' 
#' ----------------------
#' 
#' ## References
