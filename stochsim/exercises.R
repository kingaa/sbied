#' ---
#' title: "Worked solutions to basic exercises"
#' author: "Carles Bretó"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 4
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
#' ---
#' 
#' [Licensed under the Creative Commons Attribution-NonCommercial license](http://creativecommons.org/licenses/by-nc/4.0/).
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' 
## ----include=FALSE-------------------------------------------------------
source("stochsim.R")

#' 
#' ## Basic Exercise: Explore the SIR model
#' 
#' The simulated data seem to fail to capture different aspects of the data. In particular, the simulated data appear to peak substantially later (if at all) than the observed data. This in turn results in simulated valleys that arrive (again, if at all) much later. To attempt to simulate data for which the observed data is a more plausible realization, one could try increasing the force of infection.
#' 
## ------------------------------------------------------------------------
sir %>%
  simulate(params=c(Beta=2.5,gamma=1,rho=0.9,N=2600),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=day,y=B,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' Increasing the force of infection improves in one direction but now the peak is too tall. To counteract this, one could try reducing the population size.
#' 
## ------------------------------------------------------------------------
sir %>%
  simulate(params=c(Beta=2.5,gamma=1,rho=0.9,N=1500),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=day,y=B,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' And even perhaps increase the duration of immunity.
#' 
## ------------------------------------------------------------------------
sir %>% simulate(params=c(Beta=2.5,gamma=1.5,rho=0.9,N=1500),
  nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=day,y=B,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' ## Basic Exercise: The SEIR model
#' 
#' The existing code may be modified as follows:
#' 
## ------------------------------------------------------------------------
seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

seir_init <- Csnippet("
  S = N-1;
  E = 1;
  I = 0;
  R = 0;
  H = 0;
")

bsflu %>%
  pomp(times="day",t0=0,
    rprocess=euler(seir_step,delta.t=1/6),
    rinit=seir_init,
    paramnames=c("N","Beta","mu_EI","gamma","rho"),
    statenames=c("S","E","I","R","H"),
    accumvars="H",
    dmeasure=Csnippet("lik = dbinom(B,H,rho,give_log);"),
    rmeasure=Csnippet("B = rbinom(H,rho);")
  ) -> seir

#' 
#' One possibility is to split the original rate $\mu_{SI}$ into $\mu_{SE}$ and $\mu_{EI}=\gamma$. 
#' 
## ------------------------------------------------------------------------
seir %>%
  simulate(params=c(Beta=2.5,mu_EI=0.75,gamma=0.75,rho=0.9,N=1500),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=day,y=B,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' Again one can increase the force of infection. 
#' 
## ------------------------------------------------------------------------
seir %>% simulate(params=c(Beta=15,mu_EI=0.75,gamma=0.75,rho=0.9,N=1500),
  nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=day,y=B,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' 
#' ----------------------------
#' 
#' ## [Back to Stochastic Simulation lesson](./stochsim.html)
#' ## [**R** codes for this document](http://raw.githubusercontent.com/kingaa/sbied/master/stochsim/exercises.R)
#' 
#' ----------------------------
#' 
#' ## References
