#' ---
#' title: "SEIR3 example"
#' author: "Aaron A. King"
#' output:
#'   html_document:
#'     toc: yes
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
## ----prelims,include=FALSE,purl=TRUE,cache=FALSE-------------------------
library(pomp)
stopifnot(packageVersion("pomp")>="1.6")
library(ggplot2)
theme_set(theme_bw())
library(plyr)
library(reshape2)
library(magrittr)
options(stringsAsFactors=FALSE)
set.seed(1221234211)

#' 
#' ----------------------------
#' 
#' ## Model formulation
#' 
#' Formulate a model with a latent class and both confinement and convalescent stages.
#' Implement it in **pomp** using a compartmental model like that diagrammed below.
#' You will have to give some thought to just how to model the relationship between the data ($B$ and $C$) and the state variables.
#' 
## ----seir3_model---------------------------------------------------------
read.table("http://kingaa.github.io/sbied/stochsim/bsflu_data.txt") -> bsflu

rproc <- Csnippet("
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
  R2 += t4 - t5;
")

init <- Csnippet("
  S  = 762;
  E  = 0;
  I  = 1;
  R1 = 0;
  R2 = 0;
")

dmeas <- Csnippet("
  lik = dpois(B,rho*R1+1e-6,give_log);
")

rmeas <- Csnippet("
  B = rpois(rho*R1+1e-6);
")

bsflu %>%
    subset(select=-C) %>%
    pomp(times="day",t0=-6,
         rprocess=euler.sim(rproc,delta.t=1/5),
         initializer=init,rmeasure=rmeas,dmeasure=dmeas,
         statenames=c("S","E","I","R1","R2"),
         paramnames=c("Beta","mu_E","mu_I","mu_R1","mu_R2","rho")
         ) -> flu

#' 
#' How many parameters can reasonably be fixed?
#' How many must be estimated?
#' Obtain some ballpark estimates of the parameters and simulate to see if you can plausibly explain the data as a realization of this model.
#' 
#' 
## ----simulations---------------------------------------------------------
coef(flu) <- c(Beta=6,mu_E=0.5,mu_I=2,mu_R1=0.2,mu_R2=0.5,rho=0.9)

flu %>%
    simulate(nsim=20,as.data.frame=TRUE,include.data=TRUE) %>%
    subset(select=c(time,B,sim)) %>%
    ggplot(aes(x=time,y=B,color=(sim=="data"),group=sim))+
    geom_line()+
    guides(color=FALSE)


#' 
#' ## Using the particle filter
#' 
## ----pfilter1------------------------------------------------------------
flu %>% pfilter(Np=1000) -> pf
logLik(pf)


#' 
## ----pfilter2------------------------------------------------------------
library(foreach)
library(doParallel)

registerDoParallel()

bake(file="pfilter2.rds",seed=594717807L,
     kind="L'Ecuyer-CMRG",
     {
         foreach (nfilt=c(10,100,1000), .combine=rbind,
                  .options.multicore=list(set.seed=TRUE)) %:%
         foreach (Np=c(1000,10000,100000), .combine=rbind) %:%
         foreach (i=1:nfilt, .combine=rbind) %dopar% {
             flu %>% pfilter(Np=Np) %>% logLik() -> ll
             data.frame(nfilt=nfilt,Np=Np,loglik=ll)
         }
     }
     ) -> lls

registerDoSEQ()

lls %>%
  ggplot(aes(x=Np,y=loglik,fill=ordered(nfilt),group=interaction(nfilt,Np)))+
  geom_violin(draw_quantiles=0.5)+
  scale_x_log10(breaks=unique(lls$Np))+
  labs(fill="nfilt")

#' 
## ------------------------------------------------------------------------
fixed_params <- with(bsflu,c(mu_R1=1/(sum(B)/512),mu_R2=1/(sum(C)/512)))

#' 
#' ----------------------------
#' 
#' ## [Back to course homepage](http://kingaa.github.io/sbied)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/seir3/seir3.R)
#' 
#' ----------------------------
#' 
#' ## References
