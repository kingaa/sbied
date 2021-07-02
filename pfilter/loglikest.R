params <-
list(prefix = "llest/")

library(plyr)
library(tidyverse)
library(pomp)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="3.0")
set.seed(1221234211)

library(tidyverse)
library(pomp)

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,rho*H,give_log);
")

rmeas <- Csnippet("
  reports = rnbinom_mu(k,rho*H);
")

read_csv("https://kingaa.github.io/sbied/pfilter/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) %>%
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","mu_IR","eta","rho","k","N"),
    params=c(Beta=15,mu_IR=0.5,rho=0.5,k=10,eta=0.06,N=38000)
  ) -> measSIR

measSIR %>%
  simulate(nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,color=(.id=="data"),group=.id))+
  geom_line()+
  guides(color=FALSE)

measSIR %>% pfilter(Np=1000) -> pf
logLik(pf)

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
##     measSIR %>% pfilter(Np=Np) -> pf
##     logLik(pf) -> ll
##     data.frame(nfilt=nfilt,Np=Np,loglik=ll)
##   } -> lls
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
        measSIR %>% pfilter(Np=Np) -> pf
        logLik(pf) -> ll
        data.frame(nfilt=nfilt,Np=Np,loglik=ll)
      } -> lls
    registerDoSEQ()
    lls
 }) -> lls

lls %>%
  ggplot(aes(x=Np,y=loglik,fill=ordered(nfilt),
    group=interaction(nfilt,Np)))+
  geom_violin(draw_quantiles=c(0.1,0.5,0.9),alpha=0.7)+
  scale_x_log10(breaks=unique(lls$Np))+
  labs(fill="nfilt")
