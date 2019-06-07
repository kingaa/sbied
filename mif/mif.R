#' ---
#' title: "Iterated filtering: principles and practice"
#' author: "Edward Ionides and Aaron A. King"
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
#' Please share and remix non-commercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' Produced with **R** version `r getRversion()` and **pomp** version `r packageVersion("pomp")`.
#' 

## ----prelims,include=FALSE,purl=TRUE,cache=FALSE-------------------------
library(tidyverse)
library(pomp)
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="2.1")
theme_set(theme_bw())
set.seed(557976883)

#' 
#' -----------
#' 
#' ----------
#' 
#' ## Introduction
#' 
#' - This tutorial covers likelihood estimation via the method of iterated filtering.
#' 
#' - It presupposes familiarity with building partially observed Markov process (POMP) objects in the **R** package **pomp** [@King2016]. 
#' 
#' - This tutorial follows on from the [topic of carrying out particle filtering (also known as sequential Monte Carlo) via `pfilter` in **pomp**](../pfilter/pfilter.html). 
#' 
#' <br>
#' 
#' -----
#' 
#' ----
#' 
#' ## Objectives
#' 
#' 1. Review the available options for inference on POMP models, to put iterated filtering in context.
#' 
#' 1. Understand how iterated filtering algorithms carry out repeated particle filtering operations, with randomly perturbed parameter values, in order to maximize the likelihood.
#' 
#' 1. Gain experience carrying out statistical investigations using iterated filtering in a relatively simple situation (fitting an SIR model to data from a measles outbreak).
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ## Classification of statistical methods for POMP models
#' 
#' - Many, many statistical methods have been proposed for inference on POMP models [@He2010;@King2016].
#' 
#' - The volume of research indicates both the importance and the difficulty of the problem.
#' 
#' - Let's start by considering three criteria to categorize inference methods:
#'     + the plug-and-play property.
#'     + full-information or feature-based
#'     + frequentist or Bayesian.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Plug-and-play (also called simulation-based) methods
#' 
#' - Inference methodology that calls `rprocess` but not `dprocess` is said to be *plug-and-play*.
#'   All popular modern Monte Carlo methods fall into this category. 
#' - Simulation-based is equivalent to plug-and-play. 
#'     + Historically, simulation-based meant simulating forward from initial conditions to the end of the time series. 
#'     + However, particle filtering methods instead consider each observation interval sequentially.
#'       They carry out multiple, carefully selected, simulations over each interval.
#'     + Plug-and-play methods can call `dmeasure`.
#'       A method that uses only `rprocess` and `rmeasure` is called "doubly plug-and-play".
#' - Two *non-plug-and-play* methods (expectation-maximization (EM) algorithms and Markov chain Monte Carlo (MCMC)) have theoretical convergence problems for nonlinear POMP models.
#'   The failures of these two workhorses of statistical computation have prompted development of alternative methodology.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Full-information and feature-based methods
#' 
#' - *Full-information* methods are defined to be those based on the likelihood function for the full data (i.e., likelihood-based frequentist inference and Bayesian inference).
#' - *Feature-based* methods either consider a summary statistic (a function of the data) or work with an an alternative to the likelihood.
#' - Asymptotically, full-information methods are statistically efficient and feature-based methods are not.
#' - In some cases, loss of statistical efficiency might be an acceptable tradeoff for advantages in computational efficiency.
#' - However:
#' 	+ Good low-dimensional summary statistics can be hard to find. 
#' 	+ When using statistically inefficient methods, it can be hard to know how much information you are losing. 
#' 	+ Intuition and scientific reasoning can be inadequate tools to derive informative low-dimensional summary statistics [@shrestha11;@ionides11-statSci].
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Bayesian and frequentist methods
#' 
#' - Recently, plug-and-play Bayesian methods have been discovered:
#'     + particle Markov chain Monte Carlo (PMCMC) [@Andrieu2010].
#'     + approximate Bayesian computation (ABC) [@Toni2009].
#' - Prior belief specification is both the strength and weakness of Bayesian methodology:
#'     + The likelihood surface for nonlinear POMP models often contains nonlinear ridges and variations in curvature. 
#'     + These situations bring into question the appropriateness of independent priors derived from expert opinion on marginal distributions of parameters.
#'     + They also are problematic for specification of "flat" or "uninformative" prior beliefs.
#'     + Expert opinion can be treated as data for non-Bayesian analysis.
#'       However, our primary task is to identify the information in the data under investigation, so it can be helpful to use methods that do not force us to make our conclusions dependent on quantification of prior beliefs.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Full-information, plug-and-play, frequentist methods
#' 
#' - Iterated filtering methods [@ionides06;@ionides15] are the only currently available, full-information, plug-and-play, frequentist methods for POMP models.
#' - Iterated filtering methods have been shown to solve likelihood-based inference problems for epidemiological situations which are computationally intractable for available Bayesian methodology [@ionides15].
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Table of POMP inference methodologies
#' 
#' <style>
#' td, th {
#' 	vertical-align: center;
#' 	horizontal-align: center;
#' 	text-align: left;
#' 	padding: 10px;
#' 	border: 1px solid black;
#' }
#' .rowhead {
#' 	transform: rotate(-90deg);
#' 	text-align: center;
#' 	padding-top: 10px;
#' 	padding-bottom: 10px;
#' }
#' </style>
#' 
#' <table>
#' <tr><th></th><th></th><th>Frequentist</th><th>Bayesian</th></tr>
#' <tr><th rowspan="4" class="rowhead">Plug-and-play</th><th>Full-information</th><td>iterated filtering</td><td>particle MCMC</td></tr>
#' <tr><th rowspan="3">Feature-based</th><td>simulated moments</td><td>ABC</td></tr>
#' <tr><td>synthetic likelihood (SL)</td><td>SL-based MCMC</td></tr>
#' <tr><td>nonlinear forecasting</td><td>&nbsp;</td></tr>
#' <tr><th rowspan="4" class="rowhead">Not plug-and-play</th><th rowspan="2">Full-information</th><td>EM algorithm</td><td>MCMC</td></tr>
#' <tr><td>Kalman filter</td><td>&nbsp;</td></tr>
#' <tr><th rowspan="2">Feature-based</th><td>Yule-Walker<sup>1</sup></td><td>extended Kalman filter<sup>2</sup></td></tr>
#' <tr><td>extended Kalman filter<sup>2</sup></td><td>&nbsp;</td></tr>
#' </table>
#' 
#' 1. Yule-Walker is the method of moments for ARMA, a linear Gaussian POMP.
#' 1. The Kalman filter gives the exact likelihood for a linear Gaussian POMP.
#'    The extended Kalman filter gives an approximation for nonlinear models that can be used for quasi-likelihood or quasi-Bayesian inference.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ## An iterated filtering algorithm (IF2)
#' 
#' We focus on the IF2 algorithm of @ionides15.
#' In this algorithm:
#' 
#' - Each iteration consists of a particle filter, carried out with the parameter vector, for each particle, doing a random walk.
#' - At the end of the time series, the collection of parameter vectors is recycled as starting parameters for the next iteration.
#' - The random-walk variance decreases at each iteration.
#' 
#' In theory, this procedure converges toward the region of parameter space maximizing the maximum likelihood.  
#' In practice, we can test this claim on examples.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### IF2 algorithm pseudocode
#' 
#' __Input__:  
#' Simulators for $f_{X_0}(x_0;\theta)$ and $f_{X_n|X_{n-1}}(x_n| x_{n-1}; \theta)$;  
#' evaluator for $f_{Y_n|X_n}(y_n| x_n;\theta)$;  
#' data, $y^*_{1:N}$ 
#' 
#' __Algorithmic parameters__:  
#' Number of iterations, $M$;  
#' number of particles, $J$;  
#' initial parameter swarm, $\{\Theta^0_j, j=1,\dots,J\}$;  
#' perturbation density, $h_n(\theta|\varphi;\sigma)$;  
#' perturbation scale, $\sigma_{1{:}M}$ 
#' 
#' __Output__:  
#' Final parameter swarm, $\{\Theta^M_j, j=1,\dots,J\}$ 
#' 
#' __Procedure__:
#' 
#' 1. $\quad$ For $m$ in $1{:} M$
#' 2. $\quad\quad\quad$ $\Theta^{F,m}_{0,j}\sim h_0(\theta|\Theta^{m-1}_{j}; \sigma_m)$ for $j$ in $1{:} J$
#' 3. $\quad\quad\quad$ $X_{0,j}^{F,m}\sim f_{X_0}(x_0 ; \Theta^{F,m}_{0,j})$ for $j$ in $1{:} J$
#' 4. $\quad\quad\quad$ For $n$ in $1{:} N$
#' 5. $\quad\quad\quad\quad\quad$ $\Theta^{P,m}_{n,j}\sim h_n(\theta|\Theta^{F,m}_{n-1,j},\sigma_m)$ for $j$ in $1{:} J$
#' 6. $\quad\quad\quad\quad\quad$ $X_{n,j}^{P,m}\sim f_{X_n|X_{n-1}}(x_n | X^{F,m}_{n-1,j}; \Theta^{P,m}_{n,j})$ for $j$ in $1{:} J$
#' 7. $\quad\quad\quad\quad\quad$ $w_{n,j}^m = f_{Y_n|X_n}(y^*_n| X_{n,j}^{P,m} ; \Theta^{P,m}_{n,j})$ for $j$ in $1{:} J$
#' 8. $\quad\quad\quad\quad\quad$ Draw $k_{1{:}J}$ with $P[k_j=i]=  w_{n,i}^m\Big/\sum_{u=1}^J w_{n,u}^m$
#' 9.  $\quad\quad\quad\quad\quad$ $\Theta^{F,m}_{n,j}=\Theta^{P,m}_{n,k_j}$ and $X^{F,m}_{n,j}=X^{P,m}_{n,k_j}$ for $j$ in $1{:} J$
#' 10. $\quad\quad\quad$ End For
#' 11. $\quad\quad\quad$ Set $\Theta^{m}_{j}=\Theta^{F,m}_{N,j}$ for $j$ in $1{:} J$
#' 12. $\quad$ End For
#' 
#' __Remarks__:
#' 
#' - The $N$ loop (lines 4 through 10) is a basic particle filter applied to a model with stochastic perturbations to the parameters.
#' - The $M$ loop repeats this particle filter with decreasing perturbations.
#' - The superscript $F$ in $\Theta^{F,m}_{n,j}$ and $X^{F,m}_{n,j}$ denote solutions to the _filtering_ problem, with the particles $j=1,\dots,J$ providing a Monte Carlo representation of the conditional distribution at time $n$ given data $y^*_{1:n}$ for filtering iteration $m$.
#' - The superscript $P$ in $\Theta^{P,m}_{n,j}$ and $X^{P,m}_{n,j}$ denote solutions to the _prediction_ problem, with the particles $j=1,\dots,J$ providing a Monte Carlo representation of the conditional distribution at time $n$ given data $y^*_{1:n-1}$ for filtering iteration $m$.
#' - The _weight_ $w^m_{n,j}$ gives the likelihood of the data at time $n$ for particle $j$ in filtering iteration $m$.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Analogy with evolution by natural selection
#' 
#' - The parameters characterize the *genotype*.
#' - The swarm of particles is a *population*.
#' - The likelihood, a measure of the compatibility between the parameters and the data, is the analogue of *fitness*.
#' - Each successive observation is a new *generation*.
#' - Since particles reproduce in each generation in proportion to their likelihood, the particle filter acts like *natural selection*.
#' - The artificial perturbations augment the "genetic" variance and therefore correspond to *mutation*.
#' - IF2 increases the *fitness* of the population of particles.
#' - However, because our scientific interest focuses on the model without the artificial perturbations, we decrease the intensity of the latter with successive iterations.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ## Applying IF2 to the Consett measles outbreak
#' 
#' To apply IF2 to a relatively simple epidemiological example, we consider fitting a stochastic SIR model to data from a measles outbreak in the small English town of Consett.
#' Reports consist of the number of cases for each week of the year.
#' The population of the town was approximately 38,000 and, since the outbreak is confined to less than one year, we will ignore births and deaths.
#' The data are available on the course website:
#' 
## ----load_data-----------------------------------------------------------
library(tidyverse)
read_csv("https://kingaa.github.io/sbied/mif/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) %>%
  filter(week<=42) -> dat

dat %>%
  ggplot(aes(x=week,y=reports))+
  geom_line()

#' 
#' Our model is a variation on a basic SIR Markov chain, with state $X(t)=(S(t),I(t),R(t))$ giving the numbers of individuals in the susceptible and infectious categories, and three stages of recovery.
#' 

#' 
#' We'll assume that there was a single infection in week 0, i.e., that $I(0)=1$.
#' Our Markov transmission model is that each individual in $S$ transitions to $I$ at rate $\lambda=\beta\,I(t)/N$;
#' and each individual in $I$ transitions at rate $\gamma$ to $R$.
#' Therefore, $1/\gamma$ is the mean infectious period.
#' All rates will have units wk^-1^. 
#' 
#' This model has limitations and weaknesses, but writing down and fitting a model is a starting point for data analysis, not an end point.
#' In particular, having fit one model, one should certainly examine alternative models.
#' For example, one could include a latency period for infections, or one could modify the model to give a better description of the diagnosis, bed-confinement, and convalescence processes.
#' 
#' Notice that we do not need to track $R$ since this variable has consequences neither for the dynamics of the state process nor for the data.
#' 
#' Now, as before, we code the model using C snippets:
#' 
## ----csnippets-----------------------------------------------------------
library(tidyverse)
library(pomp)

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
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
  lik = dbinom(reports,H,rho,give_log);
")

rmeas <- Csnippet("
  reports = rbinom(H,rho);
")

dat %>%
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/6),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    statenames=c("S","I","H"),
    paramnames=c("Beta","gamma","eta","rho","N")
  ) -> measSIR

#' 
#' <!---
#' Note that, in our measurement model, we've added a small positive number ($10^{-6}$) to the expected number of cases.
#' Why is this useful?
#' What complications does it introduce in the interpretation of results?
#' --->
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
params <- c(Beta=15,gamma=0.5,rho=0.5,eta=0.06,N=38000)

#' 
#' Now to run and plot some simulations:
## ----init_sim------------------------------------------------------------
measSIR %>%
  simulate(params=params,nsim=10,format="data.frame") -> y


#' 
#' Before engaging in iterated filtering, it is a good idea to check that the basic particle filter is working since we can't iterate something unless we can run it once!
#' The simulations above check the `rprocess` and `rmeasure` codes;
#' the particle filter depends on the `rprocess` and `dmeasure` codes and so is a check of the latter.
#' 
## ----init_pfilter--------------------------------------------------------
measSIR %>%
  pfilter(Np=1000,params=params) -> pf


#' 
#' The above plot shows the data (`reports`), along with the *effective sample size* of the particle filter (`ess`) and the log likelihood of each observation conditional on the preceding ones (`cond.logLik`).
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Setting up the estimation problem.
#' 
#' Let's assume that the population size, $N$, is known accurately.
#' We'll fix that parameter.
#' 
## ----fixed_params--------------------------------------------------------
fixed_params <- params[c("N")]; fixed_params

#' 
#' We will estimate $\beta$, $\gamma$, $\eta$, and $\rho$.
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
#' We proceed to carry out replicated particle filters at an initial guess of $\beta=`r params["Beta"]`$, $\gamma=`r params["gamma"]`$, $\eta=`r params["eta"]`$, and $\rho=`r params["rho"]`$.
#' 
## ----pf------------------------------------------------------------------
library(doRNG)
registerDoRNG(625904618)

tic <- Sys.time()

foreach(i=1:10,.packages='pomp') %dopar% {
  measSIR %>% pfilter(params=params,Np=10000)
} -> pf

(L_pf <- logmeanexp(sapply(pf,logLik),se=TRUE))

toc <- Sys.time()

#' 
#' In `r round(toc-tic,2)` seconds, using `r min(getDoParWorkers(),length(pf))` cores, we obtain an unbiased likelihood estimate of `r round(L_pf[1],1)` with a Monte Carlo standard error of `r signif(L_pf[2],2)`.
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Building up a picture of the likelihood surface
#' 
#' Given a model and a set of data, the likelihood surface is well defined, though it may be difficult to visualize.
#' We can develop a progressively more complete picture of this surface by storing likelihood estimates whenever we compute them.
#' In particular, it is a very good idea to set up a database within which to store the likelihood of every point for which we have an estimated likelihood.
#' This will become larger and more complete as our parameter-space search goes on and will be a basis for a variety of explorations.
#' At this point, we've computed the likelihood at a single point.
#' Let's store this point, together with the estimated likelihood and our estimate of the standard error on that likelihood, in a CSV file:
## ----init_csv------------------------------------------------------------
pf[[1]] %>% coef() %>% bind_rows() %>%
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) %>%
  write_csv(path="measles_params.csv")

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
#' Since $\beta$ and $\gamma$ will be estimated on the log scale, and we expect that multiplicative perturbations of these parameters will have roughly similar effects on the likelihood, we'll use a perturbation size of $0.02$, which we imagine will have a small but non-negligible effect.
#' For simplicity, we'll use the same perturbation size on $\rho$.
#' We fix `cooling.fraction.50=0.5`, so that after 50 `mif2` iterations, the perturbations are reduced to half their original magnitudes.
#' 
## ----local_search--------------------------------------------------------
registerDoRNG(482947940)
bake(file="local_search.rds",{
  foreach(i=1:20,.packages='pomp',.combine=c) %dopar%  
  {
    measSIR %>%
    mif2(
      params=params,
      partrans=parameter_trans(log=c("Beta","gamma"),logit=c("rho","eta")),
      paramnames=c("Beta","gamma","rho","eta"),
      Np=2000,
      Nmif=50,
      cooling.fraction.50=0.5,
      rw.sd=rw.sd(Beta=0.02,gamma=0.02,rho=0.02,eta=0.02)
    )
  }
}) -> mifs_local

#' 
#' We obtain some diagnostic plots with the `plot` command applied to `mifs_local`.
#' Here is a way to get a prettier version:
#' 

#' 
#' No filtering failures (`nfail`) are generated after about 15 iterations, which is comforting.
#' In general, we expect to see filtering failures whenever our initial guess (`start`) is incompatible with one or more of the observations.
#' Filtering failures at the MLE are an indication that the model, at its best, is incompatible with one or more of the data.
#' 
#' We see that the likelihood eventually increases as the iterations proceed, though there is considerable variability due to the poorness of our starting guess and the stochastic nature of this Monte Carlo algorithm.
#' 
#' Although the filtering carried out by `mif2` in the final filtering iteration generates an approximation to the likelihood at the resulting point estimate, this is not usually good enough for reliable inference.
#' Partly, this is because parameter perturbations are applied in the last filtering iteration, so that the likelihood shown here is not identical to that of the model of interest.
#' Partly, this is because `mif2` is usually carried out with fewer particles than are needed for a good likelihood evaluation:
#' the errors in `mif2` average out over many iterations of the filtering.
#' Therefore, we evaluate the likelihood, together with a standard error, using replicated particle filters at each point estimate:
#' 
## ----lik_local-----------------------------------------------------------
registerDoRNG(900242057)
bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.packages='pomp',.combine=rbind) %dopar% 
  {
    evals <- replicate(10, logLik(pfilter(mf,Np=20000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
        bind_cols(loglik=ll[1],loglik.se=ll[2])
  }
}) -> results


#' 
#' This investigation took `r round(attr(mifs_local,"system.time")["elapsed"],0)` sec for the maximization and `r round(t_local["elapsed"],0)` sec for the likelihood evaluation.
#' These repeated stochastic maximizations can also show us the geometry of the likelihood surface in a neighborhood of this point estimate:
#' 

#' 
#' This plot shows hints of ridges in the likelihood surface (cf. the $\beta$-$\eta$ and $\beta$-$\rho$ panels).
#' However, the sampling is still too sparse to give a clear picture.
#' 
#' We add these newly explored points to our database:
## ----local_database------------------------------------------------------
results %>%
  write_csv(path="measles_params.csv",append=TRUE)

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
#' $\beta\in (10,80)$, $\gamma\in (0.3,2)$, $\rho\in (0.2,0.9)$, $\eta\in (0,0.2)$.
#' 
#' We are now ready to carry out likelihood maximizations from diverse starting points.
#' 
## ----global_search-------------------------------------------------------
registerDoRNG(1270401374)
guesses <- runifDesign(
  lower=c(Beta=10,gamma=0.3,rho=0.2,eta=0),
  upper=c(Beta=80,gamma=2,rho=0.9,eta=0.2),
  nseq=300
)

mf1 <- mifs_local[[1]]

bake(file="global_search.rds",{
  foreach(guess=iter(guesses,"row"), 
    .packages='pomp', 
    .combine=rbind,
    .export=c("mf1","fixed_params")
  ) %dopar% 
  {
    mf1 %>%
      mif2(params=c(unlist(guess),fixed_params)) %>%
      mif2(Nmif=100) -> mf
    ll <- replicate(10,mf %>% pfilter(Np=100000) %>% logLik())
    ll <- logmeanexp(ll,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
        bind_cols(loglik=ll[1],loglik.se=ll[2])
  }
}) -> results


## ------------------------------------------------------------------------
results %>%
  write_csv(path="measles_params.csv",append=TRUE)

#' 
#' The above codes run one search from each of `r nrow(guesses)` starting values.
#' Each search consists of an initial run of `r nrow(traces(mf1))` IF2 iterations, followed by another 100 iterations.
#' These codes exhibit a general **pomp** behavior:
#' re-running a command on an object (i.e., `mif2` on `mf1`) created by the same command preserves the algorithmic arguments.
#' In particular, running `mif2` on the result of a `mif2` computation re-runs IF2 from the endpoint of the first run.
#' In the second computation, by default, all algorithmic parameters are preserved;
#' here we overrode the default choice of `Nmif`.
#' 
#' Following the `mif2` computations, the particle filter is used to evaluate the likelihood, as before.
#' In contract to the local-search codes above, here we return only the endpoint of the search, together with the likelihood estimate and its standard error in a named vector.
#' The best result of this search had a likelihood of `r round(max(results$loglik),1)` with a standard error of `r round(results$loglik.se[which.max(results$loglik)],2)`.
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
#' 
## ----eta_profile---------------------------------------------------------
registerDoRNG(830007657)

read_csv("measles_params.csv") %>%
  filter(loglik>max(loglik)-20,loglik.se<2) %>%
  sapply(range) -> box
box

guesses <- profileDesign(
  eta=seq(0.01,0.85,length=20),
  lower=box[1,c("Beta","gamma","rho")],
  upper=box[2,c("Beta","gamma","rho")],
  nprof=15
)

mf1 <- mifs_local[[1]]

bake(file="eta_profile.rds",{
  foreach(guess=iter(guesses,"row"), 
    .packages='pomp', 
    .combine=rbind,
    .export=c("mf1","fixed_params")
  ) %dopar% 
  {
    mf1 %>%
      mif2(params=c(unlist(guess),fixed_params),
           rw.sd=rw.sd(Beta=0.02,gamma=0.02,rho=0.02)) %>%
      mif2(Nmif=100,cooling.fraction.50=0.3) -> mf
    ll <- replicate(10,mf %>% pfilter(Np=100000) %>% logLik())
    ll <- logmeanexp(ll,se=TRUE)
    mf %>% coef() %>% bind_rows() %>%
        bind_cols(loglik=ll[1],loglik.se=ll[2])
  }
}) -> results


## ------------------------------------------------------------------------
results %>%
  write_csv(path="measles_params.csv",append=TRUE)

## ------------------------------------------------------------------------
read_csv("measles_params.csv") %>%
  filter(loglik>max(loglik)-50) -> all

pairs(~loglik+Beta+gamma+eta+rho, data=all,
      col=ifelse(all$type=="guess",grey(0.5),"red"),pch=16)

all %>%
  filter(loglik>max(loglik)-10) %>%
  ggplot(aes(x=eta,y=loglik))+
  geom_point()

#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### Exercises
#' 
#' --------------------------
#' 
#' #### Basic Exercise: Fitting the SEIR^3^ model
#' 
#' Following the template above, estimate the parameters and likelihood of the SEIR^3^ model you implemented in the earlier lessons.
#' Specifically, first conduct a local search and then a global search using the multi-stage, multi-start method displayed above.
#' 
#' How does the maximized likelihood compare with what we obtained for the SIR^3^ model?
#' How do the parameter estimates differ?
#' 
#' You will need to tailor the intensity of your search to the computational resources at your disposal.
#' In particular, choose the number of starts, number of particles employed, and the number of IF2 iterations to perform in view of the size and speed of your machine.
#' 
#' 
#' --------------------------
#' 
#' #### Extra Exercise: Modify the measurement model
#' 
#' The Poisson measurement model used here may not seem easy to interpret.
#' Formulate an alternative measurement model and maximize the likelihood to compare the alternative model.
#' 
#' --------------------------
#' 
#' #### Extra Exercise: Construct a profile likelihood
#' 
#' How strong is the evidence about the contact rate, $\beta$, given this model and data?
#' Use `mif2` to construct a profile likelihood.
#' Due to time constraints, you may be able to compute only a preliminary version.
#' 
#' It is also possible to profile over the basic reproduction number, $R_0=\beta /\gamma$.
#' Is this more or less well determined that $\beta$ for this model and data?
#' 
#' --------------------------
#' 
#' #### Extra Exercise: Checking the source code
#' 
#' Check the source code for the `flu` pomp object.
#' Does the code implement the model described?
#' 
#' For various reasons, it can be surprisingly hard to make sure that the written equations and the code are perfectly matched.
#' Here are some things to think about:
#' 
#' 1. Papers should be written to be readable.
#'    Code must be written to run successfully.
#'    People rarely choose to clutter papers with numerical details which they hope and believe are scientifically irrelevant.
#'    What problems can arise due to the conflict between readability and reproducibility?
#'    What solutions are available?
#' 
#' 1. Suppose that there is an error in the coding of `rprocess` and suppose that plug-and-play statistical methodology is used to infer parameters.
#'    As a conscientious researcher, you carry out a simulation study to check the soundness of your inference methodology on this model.
#'    To do this, you use `simulate` to generate realizations from the fitted model and checking that your parameter inference procedure recovers the known parameters, up to some statistical error.
#'    Will this procedure help to identify the error in `rprocess`?
#'    If not, how might you debug `rprocess`?
#'    What research practices help minimize the risk of errors in simulation code?
#' 
#' <br>
#' 
#' 
#' --------
#' 
#' -------
#' 
#' ## Choosing the algorithmic settings for IF2
#' 
#' + The initial parameter swarm, $\{ \Theta^0_j, j=1,\dots,J\}$, usually consists of $J$ identical replications of some starting parameter vector.
#' + $J$ is set to be sufficient for particle filtering.
#'   Because the addition of random perturbations acts to combat particle depletion, it is typically possible to take $J$ substantially smaller than the value needed to obtain precise likelihood estimates via `pfilter`.
#'   By the time of the last iteration ($m=M$) one should not have effective sample size close to 1. 
#' + Perturbations are usually chosen to be multivariate normal, with $\sigma_m$ being a scale factor for iteration $m$:
#' $$h_n(\theta|\varphi;\sigma) \sim N[\varphi, \sigma^2_m V_n].$$
#' + $V_n$ is usually taken to be diagonal,
#' $$ V_n = \left( \begin{array}{ccccc}
#' v_{1,n}^2 & 0 & 0 & \cdots & 0 \\
#' 0 & v_{2,n}^2 &  0 & \cdots & 0 \\
#' 0 & 0 & v_{3,n}^2 & \cdots & 0 \\
#' \vdots & \vdots & \vdots & \ddots & \vdots \\
#' 0 & 0 & 0 & \cdots & v_{p,n}^2 \end{array}\right).$$
#' + If $\theta_i$ is a parameter that affects the dynamics or observations throughout the time series, it is called a __regular parameter__, and it is often appropriate to specify $$v_{i,n} = v_i.$$
#' + If $\theta_j$ is a parameter that affects only the initial conditions of the dynamic model, it is called an __initial value parameter__ (IVP) and it is appropriate to specify $$v_{j,n} = \left\{\begin{array}{ll} v_j & \mbox{if $n=0$} \\0 & \mbox{if $n>0$} \end{array}\right.$$
#' + If $\theta_k$ is a break-point parameter that models how the system changes at time $t_q$, then $\theta_k$ is like an IVP at time $t_q$ and it is appropriate to specify $$v_{j,n} = \left\{\begin{array}{ll} v_j & \mbox{if $n=q$} \\	0 & \mbox{if $n\neq q$} \end{array}\right.$$
#' + $\sigma_{1:M}$ is called a __cooling schedule__, following a thermodynamic analogy popularized by [simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing).
#'   As $\sigma_m$ becomes small, the system cools toward a "freezing point".
#'   If the algorithm is working successfully, the freezing point should be close to the lowest-energy state of the system, i.e., the MLE.
#'   Typical choices of the cooling schedule are geometric, $\sigma_m = \alpha^m$, and hyperbolic, $\sigma_m \propto 1/(1+\alpha\,m)$.
#'   In `mif2`, the cooling schedule is parameterized by $\sigma_{50}$, the cooling fraction after 50 IF2 iterations.
#' + It is generally helpful to transform the parameters so that (on the estimation scale) they are real-valued, unconstrained, and have uncertainty on the order of 1 unit.
#'   For example, one typically takes a logarithmic transformation of positive parameters and a logistic transformation of $[0,1]$ valued parameters.
#' + On such a scale, it is surprisingly often effective to take $$v_i \sim 0.02$$ for regular parameters (RPs) and $$v_j \sim 0.1$$ for initial value parameters (IVPs).
#' + We suppose that $\sigma_1=1$, since the scale of the parameters is addressed by the matrix $V_n$.
#'   Early on in an investigation, one might take $M=100$ and $\sigma_M=0.1$.
#'   As the investigation proceeds, consideration of diagnostic plots may suggest refinements. 
#' + It is remarkable that useful general advice exists for the choice of algorithmic parameters that should in principle be model- and data-specific.
#'   Here is one possible explanation:
#'   the precision of interest is often the second significant figure and there are often order 100 observations (10 monthly observations would be too few to fit a mechanistic model;
#'   1000 would be unusual for an epidemiological system). 
#' 
#' <br>
#' 
#' ----
#' 
#' ----
#' 
#' ### More exercises
#' 
#' --------------------------
#' 
#' #### Optional Exercise: Assessing and improving algorithmic parameters
#' 
#' Develop your own heuristics to try to improve the performance of `mif2` in the previous example.
#' Specifically, for a global optimization procedure carried out using random starting values in the specified box, let
#' $\hat\Theta_{\mathrm{max}}$ be a random Monte Carlo estimate of the resulting MLE, and let $\hat\theta$ be the true (unknown) MLE.
#' We can define the maximization error in the log likelihood to be
#' $$e = \ell(\hat\theta) - E[\ell(\hat\Theta_{\mathrm{max}})].$$
#' We cannot directly evaluate $e$, since there is also Monte Carlo error in our evaluation of $\ell(\theta)$, but we can compute it up to a known precision.
#' Plan some code to estimates $e$ for a search procedure using a computational effort of $JM=2\times 10^7$, comparable to that used for each mif computation in the global search.
#' Discuss the strengths and weaknesses of this quantification of optimization success.
#' See if you can choose $J$ and $M$ subject to this constraint, together with choices of `rw.sd` and the cooling rate, `cooling.fraction.50`, to arrive at a quantifiably better procedure.
#' Computationally, you may not be readily able to run your full procedure, but you could run a quicker version of it.
#' 
#' --------------------------
#' 
#' #### Optional Exercise: Finding sharp peaks in the likelihood surface
#' 
#' Even in this small, 3 parameter example, it takes a considerable amount of computation to find the global maximum (with values of $\beta$ around 0.004) starting from uniform draws in the specified box.
#' The problem is that, on the scale on which "uniform" is defined, the peak around $\beta\approx 0.004$ is very narrow.
#' Propose and test a more favorable way to draw starting parameters for the global search, with better scale invariance properties.
#' 
#' --------------------------
#' 
#' ## [Back to course homepage](../index.html)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/mif/mif.R)
#' 
#' ----------------------
#' 
#' ## References
