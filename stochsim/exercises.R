#' ---
#' title: "Worked solutions to basic exercises"
#' author: "Aaron King and Carles Bret&oacute;"
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
#' In the simulated outbreaks, the overall incidence is much too low, and the outbreak dies out after only a few weeks. 
#' To attempt to simulate data for which the observed data is a more plausible realization, we might try increasing the force of infection.
#' 
## ------------------------------------------------------------------------
sir %>%
  simulate(params=c(Beta=20,gamma=0.5,rho=0.5,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' Taking it farther....
#' 
## ------------------------------------------------------------------------
sir %>%
  simulate(params=c(Beta=40,gamma=0.5,rho=0.5,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' While this increases the overall incidence, the epidemic is now peaking too quickly.
#' To counteract this, we might try reducing the recovery rate.
#' 
## ------------------------------------------------------------------------
sir %>%
  simulate(params=c(Beta=40,gamma=0.2,rho=0.5,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' Additionally, we might have a look at the effects of changing the initial susceptible fraction, $\eta$.
#' It's possible to get something not too awful to contemplate by just manipulating $\eta$, in fact:
#' 
## ------------------------------------------------------------------------
sir %>%
  simulate(params=c(Beta=15,gamma=0.5,rho=0.5,eta=0.06,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' 
#' 
#' ## Basic Exercise: The SEIR model
#' 
#' The existing code may be modified as follows:
#' 
## ------------------------------------------------------------------------
seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-sigma*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


seir_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

sir %>%
  pomp(
    rprocess=euler(seir_step,delta.t=1/6),
    rinit=seir_init,
    paramnames=c("N","Beta","sigma","gamma","rho","eta"),
    statenames=c("S","E","I","R","H")
  ) -> seir

#' 
#' Using the rough estimate that the latent period in measles is 8--10da, we take $\sigma\sim 0.8$wk^-1^ and $\gamma\sim 1.3$wk^-1^ (so as to have roughly the same generation time as before).
#' 
## ------------------------------------------------------------------------
seir %>%
  simulate(params=c(Beta=15,sigma=0.8,gamma=1.3,rho=0.5,eta=0.06,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

#' 
#' Again one can increase the force of infection: 
#' 
## ------------------------------------------------------------------------
seir %>% 
  simulate(params=c(Beta=40,sigma=0.8,gamma=1.3,rho=0.5,eta=0.06,N=38000),
  nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
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
