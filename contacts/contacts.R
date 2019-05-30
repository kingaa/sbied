#' ---
#' title: "Case study: panel data on dynamic variation in sexual contact rates"
#' author: "Edward Ionides"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 3
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
#' ---
#' 
#' Licensed under the Creative Commons attribution-noncommercial license, http://creativecommons.org/licenses/by-nc/3.0/.
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
## ----opts,include=FALSE,cache=FALSE--------------------------------------
library(pomp)
library(panelPomp)
library(ggplot2)
theme_set(theme_bw())

library(foreach)
library(doParallel)

options(
  keep.source=TRUE,
  encoding="UTF-8"
)


registerDoParallel()
mcopts <- list(set.seed=TRUE)

#' 
#' ----------------------
#' 
#' ## Objectives
#' 
#' 1. Show how partially observed Markov process (POMP) methods can be used to understand the outcomes of a longitudinal behavioral survey.
#' 
#' 2. More broadly, discuss the use of POMP methods for panel data analysis (also known as longitudinal data).
#' 
#' 3. Introduce the R package **panelPomp** that extends **pomp** to panel data.
#' 
#' -------------
#' 
#' ------------
#' 
#' 
#' ## Introduction
#' 
#' * Basic epidemiological models suppose equal contact rates for all individuals in a population. 
#' 
#' * Sometimes these models are extended to permit rate heterogeneity between individuals. 
#' 
#' * Rate heterogeneity within individuals, i.e., dynamic behavioral change, has rarely been considered. 
#' 
#' * There have been some indications that rate heterogeneity plays a substantial role in the HIV epidemic. 
#' 
#' * @romero-severson15 investigated whether dynamic variation in sexual contact rates are a real and measurable phenomenon, by analyzing a large cohort study of HIV-negative gay men in 3 cities [@vittinghoff99]. 
#' 
#' * @romero-severson15 found evidence for dynamic variation. In a simple model for HIV, with a fully mixing population of susceptible and infected individuals, this fitted variation can help explain the observed prevalence history in the US despite the low per-contact infectivity of HIV.
#' 
#' * Here, we consider the longitudinal data from @vittinghoff99 on total sexual contacts over four consecutive 6-month periods, for the 882 men having no missing observations. 
#' 
#' * The data are available in [contacts.csv](contacts.csv). Plotted is a sample of 15 time series in the panel:
#' 
## ----data----------------------------------------------------------------
contact_data <- read.table(file="contacts.csv",header=TRUE)
matplot(t(contact_data[1:15,1:4]),
        ylab="total sexual contacts",xlab="6-month intervals", 
        type="l",xaxp=c(1,4,3))

#' 
#' We want a model that can describe all sources of variability in the data:
#' 
#' 1. Differences between individuals
#' 
#' 2. Differences within individuals over time
#' 
#' 3. Over-dispersion: variances exceeding that of a Poisson model
#' 
#' 
#' ---------------
#' 
#' ----------------
#' 
#' ## A model for dynamic variation in sexual contact rates
#' 
#' * We use the model of @romero-severson15, with each individual having a latent rate $X_i(t)$ of making contacts of a specific type.
#' 
#' * Each data point, $y_{ij}$, is the number of reported contacts for individual $i$ between time $t_{j-1}$ and $t_j$, where $i=1,\dots,882$ and $j=1,\dots,4$.
#' 
#' * The unobserved process $\{X_i(t)\}$ is connected to the data through the expected number of contacts for individual $i$ in reporting interval $j$, which we write as
#' $$C_{ij}= \alpha^{j-1}\int_{t_{j-1}}^{t_j} X_i(t)\, dt,$$
#' where $\alpha$ is an additional secular trend that accounts for the observed decline in reported contacts.
#' 
#' * A basic stochastic model for homogeneous count data would model $y_{ij}$ as a Poisson random variable with mean and variance equal to $C_{ij}$
#' [@Keeling2009].
#' 
#' * However, the variance in the data are much higher than the mean of the data [@romero-severson12-scid].
#' 
#' * Therefore, we assume that the data are negative binomially distributed [@breto11], which is a generalization of a Poisson distribution that allows for increased variance for a fixed mean, leading to the model
#' $$y_{ij}\sim \mathrm{NegBin}\, \left(C_{ij},D_{i}\right),$$
#' with mean $C_{ij}$ and  variance $C_{ij}+C_{ij}^2/D_i$.
#' 
#' * Here, $D_i$ is called the dispersion parameter, with the Poisson model being recovered in the limit as $D_i$ becomes large.
#' 
#' * The dispersion, $D_i$, can model increased variance compared to the Poisson distribution for individual contacts, but does not result in autocorrelation between measurements on an individual over time, which is observed in the data.
#' 
#' * To model this autocorrelation, we suppose that individual $i$ has behavioral episodes within which $X_i(t)$ is constant, but the individual enters new behavioral episodes at rate $R_i$. 
#' At the start of each episode, $X_i(t)$ takes a new value drawn from a Gamma distribution with mean $\mu_X$ and variance $\sigma_X$,
#' $$X_i(t)\sim \mbox{Gamma}(\mu_X, \sigma_X).$$
#' 
#' * To complete the model, we also assume Gamma distributions for $D_i$ and $R_i$,
#' $$D_i \sim \mbox{Gamma}(\mu_D, \sigma_D),$$
#' $$R_i \sim \mbox{Gamma}(\mu_R, \sigma_R).$$
#' The parameters, $\sigma_X$, $\sigma_D$ and $\sigma_R$ control individual-level differences in behavioral parameters allowing the model to encompass a wide range of sexual contact patterns.
#' 
#' * The distinction between the effects of the rate at which new behavioral episodes begin, $R_i$, and the dispersion parameter, $D_i$, is subtle since both model within-individual variability.
#' 
#' * The signal in the data about distinct behavioral episodes could be overwhelmed by a high variance in number of reported contacts resulting from a low value of $D_i$.
#' 
#' * Whether the data are sufficient to identify both $R_i$ and $D_i$ is an empirical question.
#' 
#' ---------
#' 
#' --------
#' 
#' ### Results: Consequences of dynamic behavior in a simple epidemic model
#' 
#' * Consider an SI model for HIV where the contact rates are either (a) constant; (b) vary only between individuals; (c) vary both between and within individuals. 
#' 
#' * In each case, parameterize the model by fitting the behavioral model above, and supplying per-contact infection rates from the literature. 
#' 
#' * Though this model is too simple to draw firm scientific conclusions, it does show the importance of the issue:
#' 
#' ![fig4](contacts_fig4.jpg)
#' 
#' Fig 4 of  @romero-severson15. The median of 500 simulations are shown as lines and the $75^{th}$ and $25^{th}$ quantiles are shown as gray envelopes for three parameterizations. 
#' 
#' * 'Homogeneous' (dashed line): the epidemic was simulated where $\mu_X$ is estimated by the sample mean (1.53 $\mathrm{month}^{-1}$) without any sources of between-individual or within-individual heterogeneity.
#' 
#' * 'Between Heterogeneity' (dotted line): the epidemic was simulated where $\mu_X$ is estimated by the sample mean (1.53 $\mathrm{month}^{-1}$) and $\sigma_X$ is estimated by the sample standard deviation (3.28 $\mathrm{month}^{-1}$)
#' 
#' * 'Within+Between Heterogeneity' (solid line): the epidemic was simulated where each parameter is set to the estimated maximum likelihood estimate for total contacts.
#' 
#' * For all situations, the per contact probability of transmission was set to 1/120, the average length of infection was set to 10 years, and the infection-free equilibrium population size was set to 3000. The per contact probability was selected such that the basic reproduction number in the the 'Homogeneous' case was 1.53. In the 'Homogeneous', 'Between Heterogeneity', `Within+Between Heterogeneity' cases respectively 239/500 and 172/500, 95/500 simulations died out before the 100 year mark.
#' 
#' 
#' ---------------
#' 
#' ---------------
#' 
#' ## PanelPOMP models as an extension of POMP models
#' 
#' * A PanelPOMP model consists of independent POMP models for a collection of **units**.
#' 
#' * The POMP models are tied together by shared parameters. 
#' 
#' * Here, the units are individuals in the longitudinal survey.
#' 
#' * In general, some parameters may be **unit-specific** (different for each individual) whereas others are **shared** (common to all individuals).
#' 
#' * Here, we only have shared parameters. The heterogeneities between individuals are modeled as **random effects** with distribution determined by these shared parameters.
#' 
#' * **pomp** methods were extended to PanelPOMP models by @breto19. 
#' 
#' * The main task of **panelPomp** beyond **pomp** is to handle the additional book-keeping necessitated by the unit structure.
#' 
#' * PanelPOMP models also motivate methodological developments to deal with large datasets and the high dimensional parameter vectors that can result from unit-specific parameters.
#' 
#' * A `panelPomp` object for the above contact data and model is provided by `pancon` in **panelPomp**.
#' 
## ----load-pancon---------------------------------------------------------
pompExample(pancon)
slotNames(pancon)
class(unitobjects(pancon)[[1]])

#' 
#' * The POMP models comprising the PanelPOMP model do not need to have the same observation times for each unit, or to have the same model structure. 
#' 
#' * However, if there are no shared parameters then there is no useful panel structure and the PanelPOMP model is equivalent to a list of POMP models.
#' 
#' * If we formally write a PanelPOMP as a POMP, we can use methods such as `mif2` for inference.
#' 
#' * A naive way to do inference for a PanelPOMP model as a POMP is to let an observation for the POMP be a vector of observations for all units in the PanelPOMP at that time. This gives a high-dimensional observation vector which is numerically intractable via particle filters.
#' 
#' * Instead, the `mif2` method for `panelPomp` objects concatenates the panel members into one long time series.
#' 
#' 
#' * `pancon` contains the MLE for total contacts reported in the top row of table 1 of @romero-severson15:
#' 
## ----coef----------------------------------------------------------------
coef(pancon)

#' 
#' * The implementation of the above model equations in `pancon` can be found in the `panelPomp` source code:
#' https://github.com/cbreto/panelPomp/blob/master/inst/examples/pancon.R
#' 
#' ---------
#' 
#' #### Question
#' 
#' * Describe what you think `pfilter(pancon,Np=100)` should do. 
#' 
#' * What might be the class of the resulting object? What slots might this object possess?
#' 
#' --------
#' 
#' 
## ----pfilter1------------------------------------------------------------
stew("pfilter1.rda",{
t1 <- system.time(
pf1 <- foreach(i=1:10,.packages=c('pomp','panelPomp'),.options.multicore=mcopts) %dopar% try(pfilter(pancon,Np=1000))
)
},seed=1943810296,kind="L'Ecuyer")

(loglik1 <- sapply(pf1,logLik))

#' 
#' 
#' -----------
#' 
#' ----------
#' 
#' ### Likelihood maximization using the PIF algorithm
#' 
#' * The panel iterated filtering (PIF) algorithm of @breto19 is an extension of the IF2 algorithm of @ionides15.
#' 
#' * PIF is implemented in **panelPomp** as the `mif2` method for class `panelPomp`. 
#' 
#' * Comparing `?panelPomp::mif2` with `?pomp::mif2` reveals that the only difference in the arguments
#' is that the `start` argument for `pomp::mif2` becomes `shared.start` and `specific.start` for `panelPomp::mif2`.
#' 
#' * As an example of an iterated filtering investigation, let's carry out a local search, starting at the current estimate of the MLE.
## ----mif1----------------------------------------------------------------
stew("mif1.rda",{
  t2 <- system.time(
    m2 <- foreach(i=1:10,.packages=c('pomp','panelPomp'),.options.multicore=mcopts) %dopar% try( 
      mif2(pancon,
        Nmif=50,
        Np=1000,
        cooling.fraction.50=0.1,
        cooling.type="geometric",
        transform=TRUE,
        rw.sd=rw.sd(mu_X=0.02, sigma_X=0.02, mu_D = 0.02, sigma_D=0.02,
                    mu_R=0.02, sigma_R =0.02, alpha=0.02)
      )
    )
  )
},seed=354320731,kind="L'Ecuyer")

#' 
#' * This is a relatively quick search, taking `r round(t2["elapsed"]/60,1)` minutes. 
#' 
#' * The preliminary likelihood estimated as a consequence of running `mif2` and extracted here by `sapply(m2,logLik)` does not correspond to the actual, fixed parameter, model. It is the sequential Monte Carlo estimate of the likelihood from the last filtering iteration, and therefore will have some perturbation of the parameters. Further, one typically requires fewer particles for each filtering iteration than necessary to obtain a good likelihood estimate---stochastic errors can cancel out through the filtering iterations, rather than within any one iteration. 
#' 
#' * For promising new parameter values, it is desirable to put computational effort into evaluating the likelihood sufficient to make the Monte Carlo error small compared to one log unit.
#' 
#'  
## ----mif1-lik-eval-------------------------------------------------------
stew("mif1-lik-eval.rda",{
  params_new <- coef( m2[[which.max( sapply(m2,logLik) )]] )
  t3 <- system.time( pf3 <- foreach(i=1:10,.packages=c('pomp','panelPomp'),
                 .options.multicore=mcopts) %dopar% try(
                   pfilter(pancon,shared=params_new,Np=2000))
  )
},seed=35780731,kind="L'Ecuyer")

(loglik_new <- logmeanexp(sapply(pf3,logLik),se=TRUE))

#' 
#' *  This took `r round(t3["elapsed"]/60,1)` minutes
#' 
#' * This is the direct analog of how we combined Monte Carlo replicated with POMP models.
#' 
#' * However, here we might be able to take advantage of the panel structure to do better.
#' 
#' 
#' ---------------
#' 
#' ------------------------------
#' 
#' #### Exercise: How to average repeated log likelihood estimates for a PanelPOMP model
#' 
#' 
#' * It turns out that the SMC algorithm implemented by `pfilter()` gives an unbiased Monte Carlo estimate of the likelihood.
#' 
#' * For inferential purposes, we usually work with the log likelihood. Due to Jensen's inequality, SMC has a negative bias as an estimator of the log likelihood, i.e., it systematically underestimates the log likelihood. Usually, the higher the Monte Carlo error on the likelihood, the larger this bias.
#' 
#' * Let $\hat\lambda_u^{(k)}$ be the $k$th replication of the Monte Carlo log likelihood evaluation for unit $u$.
#' 
#' * Let $\hat L_u^{(k)}=\exp\big\{\hat\lambda_u^{(k)}\big\}$ be the corresponding likelihood.
#' 
#' * Let  $\hat\lambda^{(k)}=\sum_{u=1}^U \lambda_{u}^{(k)}$ and  $\hat L^{(k)}=\exp\big\{\hat\lambda^{(k)}\big\}$ be estimates of the (log)likelihood of the entire data.
#' 
#' * Different possible estimates of the actual log likelihood $\lambda=\sum_{u=1}^U \lambda_u$ are
#' 
#' \begin{eqnarray} \hat\lambda^{[1]} &=& \frac{1}{K}\sum_{k=1}^K \hat\lambda^{(k)} 
#' \\
#'   \hat\lambda^{[2]} &=& \log \left( \frac{1}{K}\sum_{k=1}^K \hat L^{(k)} \right)  
#' \\
#'  \hat\lambda^{[3]} &=& \sum_{u=1}^U\frac{1}{K}\sum_{k=1}^K \hat\lambda^{(k)}_u  
#' \\
#'   \hat\lambda^{[4]} &=& \sum_{u=1}^U \log \left( \frac{1}{K}\sum_{k=1}^K \hat L^{(k)}_u \right)  
#' \end{eqnarray}
#' 
#' (a) Calculate empirically the mean and standard deviation of the log likelihood estimate at the MLE using 100 filter replications each with 1000 particles.
#' 
#' 
#' (b) Which of $\hat\lambda^{[1]}$,  $\hat\lambda^{[2]}$,  $\hat\lambda^{[3]}$ and  $\hat\lambda^{[4]}$ would you use, and why?
#' 
#' Some discussion of different estimators and their means and variance is given in the supplement of @breto19. 
#' 
#' 
#' --------------------------
#' 
#' ## [Back to course homepage](../index.html)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/contacts/contacts.R)
#' 
#' 
#' ----------------------
#' 
#' ## Acknowledgments
#' 
#' This tutorial is prepared for the [Simulation-based Inference for Epidemiological Dynamics](https://kingaa.github.io/sbied/) module at 11th Annual Summer Institute in Statistics and Modeling in Infectious Diseases ([SISMID 2019](https://www.biostat.washington.edu/suminst/sismid)). 
#' Previous versions were presented at SISMID in 2016, 2017, 2018.
#' 
#' ----------------------
#' 
#' ## References
