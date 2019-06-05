#' ---
#' title: 'Worked solution to "log likelihood estimation" Exercise'
#' author: "Aaron A. King"
#' output:
#'   html_document:
#'     toc: no
#'     toc_depth: 4
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

## ----prelims,purl=TRUE,cache=FALSE---------------------------------------
library(plyr)
library(tidyverse)
theme_set(theme_bw())
options(stringsAsFactors=FALSE)
library(pomp)
stopifnot(packageVersion("pomp")>"2.0.9")
set.seed(1221234211)

#' 
#' Model implementation
#' 
## ----model---------------------------------------------------------------
Csnippet("
    double N = 763;
    double t1 = rbinom(S,1-exp(-Beta*I/N*dt));
    double t2 = rbinom(E,1-exp(-mu_E*dt));
    double t3 = rbinom(I,1-exp(-mu_I*dt));
    double t4 = rbinom(R1,1-exp(-mu_R1*dt));
    double t5 = rbinom(R2,1-exp(-mu_R2*dt));
    S  -= t1;
    E  += t1 - t2;
    I  += t2 - t3;
    R1 += t3 - t4;
    R2 += t4 - t5;"
) -> rproc

Csnippet("
    S  = 762;
    E  = 0;
    I  = 1;
    R1 = 0;
    R2 = 0;"
) -> init

Csnippet("
    lik = dpois(B,rho*R1+1e-6,give_log);"
) -> dmeas

Csnippet("
    B = rpois(rho*R1+1e-6);"
) -> rmeas

bsflu %>%
  select(day,B) %>%
  pomp(times="day",t0=-6,
    rprocess=euler(rproc,delta.t=1/5),
    rinit=init,rmeasure=rmeas,dmeasure=dmeas,
    statenames=c("S","E","I","R1","R2"),
    paramnames=c("Beta","mu_E","mu_I","mu_R1","mu_R2","rho")
  ) -> flu

coef(flu) <- c(Beta=6,mu_E=0.5,mu_I=2,mu_R1=0.2,mu_R2=0.5,rho=0.9)

flu %>%
  simulate(nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=day,y=B,color=(.id=="data"),group=.id))+
  geom_line()+
  guides(color=FALSE)

#' 
#' Testing the particle filter:
#' 
## ----test----------------------------------------------------------------
flu %>% pfilter(Np=1000) -> pf
logLik(pf)

#' 
#' Now, we evaluate the dependence of log likelihood estimates 
#' on particle size and number of independent filters
#' 
## ----comps,eval=FALSE----------------------------------------------------
## library(foreach)
## library(doParallel)
## 
## registerDoParallel()
## 
## foreach (nfilt=c(10,100,1000),
##   .combine=rbind,.inorder=FALSE,
##   .options.multicore=list(set.seed=TRUE)) %:%
##   foreach (Np=c(1000,10000,100000), .combine=rbind) %:%
##   foreach (i=1:nfilt, .combine=rbind) %dopar% {
##     flu %>% pfilter(Np=Np) -> pf
##     logLik(pf) -> ll
##     data.frame(nfilt=nfilt,Np=Np,loglik=ll)
##   } -> lls

## ----comps-eval,include=FALSE--------------------------------------------
bake(file="loglikest-pfilter.rds",
  seed=594717807L,kind="L'Ecuyer-CMRG",{
    library(foreach)
    library(doParallel)
    
    registerDoParallel()
    
    foreach (nfilt=c(10,100,1000),
      .combine=rbind,.inorder=FALSE,
      .options.multicore=list(set.seed=TRUE)) %:%
      foreach (Np=c(1000,10000,100000), .combine=rbind) %:%
      foreach (i=1:nfilt, .combine=rbind) %dopar% {
        flu %>% pfilter(Np=Np) -> pf
        logLik(pf) -> ll
        data.frame(nfilt=nfilt,Np=Np,loglik=ll)
      } -> lls
    registerDoSEQ()
    lls
 }) -> lls

#' 
#' Violin plots are cute.
#' 
## ----plots---------------------------------------------------------------
lls %>%
  ggplot(aes(x=Np,y=loglik,fill=ordered(nfilt),
    group=interaction(nfilt,Np)))+
  geom_violin(draw_quantiles=c(0.1,0.5,0.9),alpha=0.7)+
  scale_x_log10(breaks=unique(lls$Np))+
  labs(fill="nfilt")

#' 
#' -----------
#' 
#' ## [Back to main lesson](./pfilter.html)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/pfilter/loglikest.R)
#' 
#' -----------
#' 
#' ## References
