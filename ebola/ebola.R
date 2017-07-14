#' ---
#' title: "Case study: forecasting Ebola"
#' author: "A. A. King, M. Domenech de Cell&egrave;s, F. M. G. Magpantay, P. Rohani, E. Ionides"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 4
#'     toc_float:
#'       collapsed: FALSE
#'       smooth_scroll: TRUE
#'     code_folding: hide
#'     highlight: haddock
#'     number_sections: FALSE
#'     df_print: kable
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
#' ---
#' 
#' \newcommand\prob[1]{\mathbb{P}\left[{#1}\right]}
#' \newcommand\expect[1]{\mathbb{E}\left[{#1}\right]}
#' \newcommand\var[1]{\mathrm{Var}\left[{#1}\right]}
#' \newcommand\dist[2]{\mathrm{#1}\left(#2\right)}
#' \newcommand\dlta[1]{{\Delta}{#1}}
#' \newcommand\scinot[2]{$#1 \times 10^{#2}$\xspace}
#' 
#' ------------------------------------
#' 
#' 
#' [Licensed under the Creative Commons Attribution-NonCommercial license](http://creativecommons.org/licenses/by-nc/4.0/).
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' Produced in **R** version `r getRversion()` using **pomp** version `r packageVersion("pomp")`.
#' 
#' ------------------------------------
#' 
#' 
## ----prelims,include=FALSE,cache=FALSE-----------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8"
  )

set.seed(594709947L)
library(ggplot2)
theme_set(theme_bw())
library(plyr)
library(reshape2)
library(magrittr)
library(pomp)
stopifnot(packageVersion("pomp")>="1.12")

#' 
#' ## Objectives
#' 
#' 1. To explore the use of POMP models in the context of an outbreak of an emerging infectious disease.
#' 1. To demonstrate the use of diagnostic probes for model criticism.
#' 1. To illustrate some forecasting methods based on POMP models.
#' 1. To provide an example that can be modified to apply similar approaches to other outbreaks of emerging infectious diseases.
#' 
#' These objectives will be achieved using a recent study [@King2015], all codes for which are available on [datadryad.org](http://dx.doi.org/10.5061/dryad.r5f30).
#' 
#' 
#' ## An emerging infectious disease outbreak
#' 
#' Let's situate ourselves at the beginning of October 2014.
#' The WHO situation report contained data on the number of cases in each of Guinea, Sierra Leone, and Liberia.
#' Key questions included:
#' 
#' 1. How fast will the outbreak unfold?
#' 1. How large will it ultimately prove?
#' 1. What interventions will be most effective?
#' 
#' As is to be expected in the case of a fast-moving outbreak of a novel pathogen in an underdeveloped country, the answers to these questions were sought in a context far from ideal:
#' 
#' - Case ascertainment is difficult and the case definition itself may be evolving.
#' - Surveillance effort is changing on the same timescale as the outbreak itself.
#' - The public health and behavioral response to the outbreak is rapidly changing.
#' 
#' The @King2015 paper focused critical attention on the economical and therefore common practice of fitting deterministic transmission models to cumulative incidence data.
#' Specifically, @King2015 showed how this practice easily leads to overconfident prediction that, worryingly, can mask their own presence.
#' The paper recommended the use of POMP models, for several reasons:
#' 
#' - Such models can accommodate a wide range of hypothetical forms.
#' - They can be readily fit to incidence data, especially during the exponential growth phase of an outbreak.
#' - Stochastic models afford a more explicit treatment of uncertainty.
#' - POMP models come with a number of diagnostic approaches built-in, which can be used to assess model misspecification.
#' 
#' ## Data and model
#' 
#' The data and **pomp** codes used to represent the transmission models are presented in [a supplement](./model.html).
#' 
#' ### Situation-report data
#' 
#' The data we focus on here are from the WHO Situation Report of 1 October 2014.
#' Supplementing these data are population estimates for the three countries..
#' 
## ----get-data,include=FALSE----------------------------------------------
base_url <- "https://kingaa.github.io/sbied/"
read.csv(paste0(base_url,"ebola/ebola_data.csv"),stringsAsFactors=FALSE,
         colClasses=c(date="Date")) -> dat
sapply(dat,class)
head(dat)

#' 
## ----popsizes,include=FALSE----------------------------------------------
populations <- c(Guinea=10628972,Liberia=4092310,SierraLeone=6190280)

#' 
## ----plot-data,echo=FALSE------------------------------------------------
dat %>%
  ggplot(aes(x=date,y=cases,group=country,color=country))+
  geom_line()

#' 
#' ### An SEIR model with gamma-distributed latent and infectious periods
#' 
#' Many of the early modeling efforts used variants on the simple SEIR model.
#' Here, we'll focus on a variant that attempts a more careful description of the duration of the latent period.
#' Specifically, this model assumes that the amount of time an infection remains latent is
#' $$\mathrm{LP} \sim \dist{Gamma}{m,\frac{1}{m\,\alpha}},$$
#' where $m$ is an integer.
#' This means that the latent period has expectation $1/\alpha$ and variance $1/(m\,\alpha)$.
#' In this document, we'll fix $m=3$.
#' 
#' We implement Gamma distributions using the so-called *linear chain trick*.
#' 
#' 
#' ![Model flow diagram.](./model_diagram.png)
#' 
## ----rproc,include=FALSE-------------------------------------------------
rSim <- Csnippet('
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
')

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

#' 
## ----skel,include=FALSE--------------------------------------------------
skel <- Csnippet('
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
')

#' 
#' The observations are modeled as a negative binomial process conditional on the number of infections.
#' That is, if $C_t$ are the reported cases at week $t$ and $H_t$ is the true incidence, then we postulate that $C_t | H_t$ is negative binomial with $$\expect{C_t|H_t} = \rho\,H_t$$ and $$\var{C_t|H_t} = \rho\,H_t\,(1+k\,\rho\,H_t).$$
#' The negative binomial process allows for overdispersion in the counts.
#' This overdispersion is controlled by parameter $k$.
#' 
## ----measmodel,include=FALSE---------------------------------------------
dObs <- Csnippet('
  double f;
  if (k > 0.0)
    f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
  else
    f = dpois(nearbyint(cases),rho*N_EI,1);
  lik = (give_log) ? f : exp(f);
')

rObs <- Csnippet('
  if (k > 0) {
    cases = rnbinom_mu(1.0/k,rho*N_EI);
  } else {
    cases = rpois(rho*N_EI);
  }')

#' 
## ----partrans,include=FALSE----------------------------------------------
toEst <- Csnippet('
  const double *IC = &S_0;
  double *TIC = &TS_0;
  TR0 = log(R0);
  Trho = logit(rho);
  Tk = log(k);
  to_log_barycentric(TIC,IC,4);
')

fromEst <- Csnippet('
  const double *IC = &S_0;
  double *TIC = &TS_0;
  TR0 = exp(R0);
  Trho = expit(rho);
  Tk = exp(k);
  from_log_barycentric(TIC,IC,4);
')

#' 
## ----pomp-construction,include=FALSE-------------------------------------
ebolaModel <- function (country=c("Guinea", "SierraLeone", "Liberia"),
                        timestep = 0.1, nstageE = 3) {

  ctry <- match.arg(country)
  pop <- unname(populations[ctry])
  nstageE <- as.integer(nstageE)

  globs <- paste0("static int nstageE = ",nstageE,";")

  dat <- subset(dat,country==ctry,select=-country)

  ## Create the pomp object
  dat %>% 
    extract(c("week","cases")) %>%
    pomp(
      times="week",
      t0=min(dat$week)-1,
      globals=globs,
      statenames=c("S",sprintf("E%1d",seq_len(nstageE)),
                   "I","R","N_EI","N_IR"),
      zeronames=c("N_EI","N_IR"),
      paramnames=c("N","R0","alpha","gamma","rho","k",
                   "S_0","E_0","I_0","R_0"),
      dmeasure=dObs, rmeasure=rObs,
      rprocess=discrete.time.sim(step.fun=rSim, delta.t=timestep),
      skeleton=vectorfield(skel),
      toEstimationScale=toEst,
      fromEstimationScale=fromEst,
      initializer=rInit) -> po
}

ebolaModel("Guinea") -> gin
ebolaModel("SierraLeone") -> sle
ebolaModel("Liberia") -> lbr

#' 
#' ## Parameter estimates
#' 
#' @King2015 estimated parameters for this model for each country.
#' A Latin hypercube design was used to initiate a large number of iterated filtering runs.
#' Profile likelihoods were computed for each country against the parameters $k$ (the measurement model overdispersion) and $R_0$ (the basic reproductive ratio).
#' Full details are given [on the datadryad.org site](http://dx.doi.org/10.5061/dryad.r5f30).
#' The results of these calculations are loaded and displayed here.
#' 
## ----load-profile,echo=FALSE---------------------------------------------
options(stringsAsFactors=FALSE)
profs <- read.csv(paste0(base_url,"/ebola/ebola-profiles.csv"))

#' 
#' The following plots the profile likelihoods.
#' The horizontal line represents the critical value of the likelihood ratio test for $p=0.01$.
## ----profiles-plots,results='hide',echo=FALSE----------------------------
library(reshape2)
library(plyr)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())

profs %>% 
  melt(id=c("profile","country","loglik")) %>%
  subset(variable==profile) %>%
  ddply(~country,mutate,dll=loglik-max(loglik)) %>%
  ddply(~country+profile+value,subset,loglik==max(loglik)) %>% 
  ggplot(mapping=aes(x=value,y=dll))+
  geom_point(color='red')+
  geom_hline(yintercept=-0.5*qchisq(p=0.99,df=1))+
  facet_grid(country~profile,scales='free')+
  labs(y=expression(l))

#' 
#' ## Diagnostics *or* Model Criticism
#' 
#' + Parameter estimation is the process of finding the parameters that are "best", in some sense, for a given model, from among the set of those that make sense for that model.
#' + Model selection, likewise, aims at identifying the "best" model, in some sense, from among a set of candidates.
#' + One can do both of these things more or less well, but no matter how carefully they are done, the best of a bad set of models is still bad.
#' 
#' + Let's investigate the model here, at its maximum-likelihood parameters, to see if we can identify problems.
#' + The guiding principle in this is that, if the model is "good", then the data are a plausible realization of that model.
#' + Therefore, we can compare the data directly against model simulations.
#' + Moreover, we can quantify the agreement between simulations and data in any way we like.
#' + Any statistic, or set of statistics, that can be applied to the data can also be applied to simulations.
#' + Shortcomings of the model should manifest themselves as discrepancies between the model-predicted distribution of such statistics and their value on the data.
#' 
#' + **pomp** provides tools to facilitate this process.
#' + Specifically, the `probe` function applies a set of user-specified *probes* or summary statistics, to the model and the data, and quantifies the degree of disagreement in several ways.
#' 
#' Let's see how this is done using the model for the Guinean outbreak.
#' 
## ----diagnostics1,echo=FALSE---------------------------------------------
library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
options(stringsAsFactors=FALSE)

profs %>%
  subset(country=="Guinea") %>%
  subset(loglik==max(loglik),
         select=-c(loglik,loglik.se,country,profile)) %>%
  unlist() -> coef(gin)

simulate(gin,nsim=20,as.data.frame=TRUE,include.data=TRUE) %>% 
  mutate(date=min(dat$date)+7*(time-1),
         is.data=ifelse(sim=="data","yes","no")) %>% 
  ggplot(aes(x=date,y=cases,group=sim,color=is.data,
         alpha=is.data))+
  geom_line()+
  guides(color=FALSE,alpha=FALSE)+
  scale_color_manual(values=c(no=gray(0.6),yes='red'))+
  scale_alpha_manual(values=c(no=0.5,yes=1))

#' 
#' The simulations appear to be growing a bit more quickly than the data.
#' Let's try to quantify this.
#' First, we'll write a function that estimates the exponential growth rate by linear regression.
#' Then, we'll apply it to the data and to 500 simulations.
#' 
#' In the following, `gin` is a `pomp` object containing the model and the data from the Guinea outbreak.
#' Below, we make use of the **magrittr** package, which provides the very useful `%>%` pipe operator.
#' A brief explanation of the syntax is [given here](https://kingaa.github.io/R_Tutorial/munging.html#the-magrittr-syntax).
#' 
## ----diagnostics-growth-rate---------------------------------------------
growth.rate <- function (y) {
  cases <- y["cases",]
  fit <- lm(log1p(cases)~seq_along(cases))
  unname(coef(fit)[2])
}
probe(gin,probes=list(r=growth.rate),nsim=500) %>% plot()

#' 
#' Do these results bear out our suspicion that the model and data differ in terms of growth rate?
#' 
#' The simulations also appear to be more highly variable around the trend than do the data.
#' 
## ----diagnostics-growth-rate-and-sd--------------------------------------
growth.rate.plus <- function (y) {
  cases <- y["cases",]
  fit <- lm(log1p(cases)~seq_along(cases))
  c(r=unname(coef(fit)[2]),sd=sd(residuals(fit)))
}
probe(gin,probes=list(growth.rate.plus),
      nsim=500) %>% plot()

#' 
#' Let's also look more carefully at the distribution of values about the trend using the 1st and 3rd quantiles.
#' Also, it looks like the data are less jagged than the simulations.
#' We can quantify this using the autocorrelation function (ACF).
#' 
## ----diagnostics2,fig.height=6-------------------------------------------
log1p.detrend <- function (y) {
  cases <- y["cases",]
  y["cases",] <- as.numeric(residuals(lm(log1p(cases)~seq_along(cases))))
  y
}

probe(gin,probes=list(
  growth.rate.plus,
  probe.quantile(var="cases",prob=c(0.25,0.75)),
  probe.acf(var="cases",lags=c(1,2,3),type="correlation",
            transform=log1p.detrend)
),nsim=500) %>% plot()

#' 
#' ### Exercise: the SEIR model for the Sierra Leone outbreak
#' 
#' Apply probes to investigate the extent to which the model is an adequate description of the data from the Sierra Leone outbreak.
#' Have a look at the probes provided with **pomp**: `?basic.probes`.
#' Try also to come up with some informative probes of your own.
#' Discuss the implications of your findings.
#' 
#' ## Forecasting using POMP models
#' 
#' To this point in the course, we've focused on using POMP models to answer scientific questions, i.e., to compare alternative hypothetical explanations for the data in hand.
#' Of course, we can also use them to make forecasts.
#' The key issues are to do with quantifying the forecast uncertainty.
#' This arises from four sources:
#' 
#' 1. measurement error
#' 1. process noise
#' 1. parametric uncertainty
#' 1. structural uncertainty
#' 
#' Here, we'll explore how we can account for the first three of these in making forecasts for the Sierra Leone outbreak.
#' 
#' We take an [*empirical Bayes*](https://en.wikipedia.org/wiki/Empirical_Bayes_method) approach:
#' 
#' - We set up a collection of parameter vectors in a neighborhood of the maximum likelihood estimate containing the region of high likelihood.
#' 
## ----forecasts1----------------------------------------------------------
library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
options(stringsAsFactors=FALSE)

set.seed(988077383L)

## forecast horizon
horizon <- 13

## Weighted quantile function
wquant <- function (x, weights, probs = c(0.025,0.5,0.975)) {
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  w <- cumsum(weights)/sum(weights)
  rval <- approx(w,x,probs,rule=1)
  rval$y
}

profs %>% 
  subset(country=="SierraLeone",
         select=-c(country,profile,loglik.se)) %>%
  subset(loglik>max(loglik)-0.5*qchisq(df=1,p=0.99)) %>%
  melt(variable.name="parameter") %>%
  ddply(~parameter,summarize,
        min=min(value),max=max(value)) %>%
  subset(parameter!="loglik") %>%
  melt(measure=c("min","max")) %>%
  acast(parameter~variable) -> ranges

params <- sobolDesign(lower=ranges[,'min'],
                      upper=ranges[,'max'],
                      nseq=20)
plot(params)

#' 
#' - We carry out a particle filter at each parameter vector, which gives us estimates of both the likelihood and the filter distribution at that parameter value.
#' - We simulate forward from each filter distribution, up to the desired forecast horizon, to give the prediction distribution for each parameter vector.
#' - We sample from these prediction distributions with probability proportional to the estimated likelihood of the parameter vector.
#' 
## ----forecasts2----------------------------------------------------------
library(foreach)
library(doParallel)
library(iterators)

registerDoParallel()

set.seed(887851050L,kind="L'Ecuyer")

foreach(p=iter(params,by='row'),
        .inorder=FALSE,
        .combine=rbind,
        .options.multicore=list(preschedule=TRUE,set.seed=TRUE)
        ) %dopar%
    {
        library(pomp)
        
        M1 <- ebolaModel("SierraLeone")

        pf <- pfilter(M1,params=unlist(p),Np=2000,save.states=TRUE)

        pf$saved.states %>%                 # latent state for each particle
            tail(1) %>%                     # last timepoint only
            melt() %>%                      # reshape and rename the state variables
            dcast(rep~variable,value.var="value") %>%
            ddply(~rep,summarize,S_0=S,E_0=E1+E2+E3,I_0=I,R_0=R) %>%
            melt(id="rep") %>%
            acast(variable~rep) -> x
        ## the final states are now stored in 'x' as initial conditions
        
        ## set up a matrix of parameters
        pp <- parmat(unlist(p),ncol(x)) 
        
        ## generate simulations over the interval for which we have data
        simulate(M1,params=pp,obs=TRUE) %>%
            melt() %>%
            mutate(time=time(M1)[time],
                   period="calibration",
                   loglik=logLik(pf)) -> calib

        ## make a new 'pomp' object for the forecast simulations
        M2 <- M1
        time(M2) <- max(time(M1))+seq_len(horizon)
        timezero(M2) <- max(time(M1))
        
        ## set the initial conditions to the final states computed above
        pp[rownames(x),] <- x
        
        ## perform forecast simulations
        simulate(M2,params=pp,obs=TRUE) %>%
            melt() %>%
            mutate(time=time(M2)[time],
                   period="projection",
                   loglik=logLik(pf)) -> proj
        
        rbind(calib,proj)
    } %>%
    subset(variable=="cases",select=-variable) %>%
    mutate(weight=exp(loglik-mean(loglik))) %>%
    arrange(time,rep) -> sims

## look at effective sample size
ess <- with(subset(sims,time==max(time)),weight/sum(weight))
ess <- 1/sum(ess^2); ess

## compute quantiles of the forecast incidence
sims %>%
    ddply(~time+period,summarize,prob=c(0.025,0.5,0.975),
          quantile=wquant(value,weights=weight,probs=prob)) %>%
    mutate(prob=mapvalues(prob,from=c(0.025,0.5,0.975),
                          to=c("lower","median","upper"))) %>%
    dcast(period+time~prob,value.var='quantile') %>%
    mutate(date=min(dat$date)+7*(time-1)) -> simq

#' 
## ----forecast-plots,echo=FALSE-------------------------------------------
simq %>% ggplot(aes(x=date))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=period),alpha=0.3,color=NA)+
  geom_line(aes(y=median,color=period))+
  geom_point(data=subset(dat,country=="SierraLeone"),
             mapping=aes(x=date,y=cases),color='black')+
  labs(y="cases")

#' 
#' --------------------------
#' 
#' ## [Back to course homepage](../index.html)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/ebola/ebola.R)
#' ## [Detailed model supplement](model.html)
#' 
#' ----------------------
#' 
#' ## References
