library(tidyverse)
library(pomp)
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="3.0")
theme_set(theme_bw())
set.seed(1350254336)





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
lik = dbinom(reports,H,rho,give_log);
")

rmeas <- Csnippet("
reports = rbinom(H,rho);
")

read_csv("https://kingaa.github.io/sbied/pfilter/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) %>%
  filter(week<=42) %>%
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","mu_IR","eta","rho","N"),
    params=c(Beta=15,mu_IR=0.5,rho=0.5,eta=0.06,N=38000)
  ) -> measSIR

measSIR %>%
  pfilter(Np=5000) -> pf
logLik(pf)

library(doParallel)
library(doRNG)
registerDoParallel()
registerDoRNG(652643293)

foreach (i=1:10, .combine=c) %dopar% {
  measSIR %>% pfilter(Np=5000)
} -> pf
logLik(pf) -> ll
logmeanexp(ll,se=TRUE)







bake(file="like-slice.rds",{
  sliceDesign(
    center=coef(measSIR),
    Beta=rep(seq(from=5,to=20,length=40),each=3),
    mu_IR=rep(seq(from=0.2,to=2,length=40),each=3)
  ) -> p
  library(doParallel)
  library(doRNG)
  
  registerDoParallel()
  registerDoRNG(108028909)
  foreach (theta=iter(p,"row"), .combine=rbind,
           .inorder=FALSE) %dopar%
    {
      library(pomp)
      measSIR %>% pfilter(params=theta,Np=5000) -> pf
      theta$loglik <- logLik(pf)
      theta
    } -> p
}) -> p

library(tidyverse)

p %>% 
  gather(variable,value,Beta,mu_IR) %>%
  filter(variable==slice) %>%
  ggplot(aes(x=value,y=loglik,color=variable))+
  geom_point()+
  facet_grid(~variable,scales="free_x")+
  guides(color=FALSE)+
  labs(x="parameter value",color="")+
  theme_bw()







bake(file="pfilter-grid1.rds",{
  expand.grid(
    Beta=rep(seq(from=10,to=30,length=40),each=3),
    mu_IR=rep(seq(from=0.4,to=1.5,length=40),each=3),
    rho=0.5,eta=0.06,N=38000
  ) -> p
  library(doParallel)
  library(doRNG)
  
  registerDoParallel()
  registerDoRNG(421776444)
  foreach (theta=iter(p,"row"), .combine=rbind,
           .inorder=FALSE) %dopar%
    {
      library(pomp)
      measSIR %>% pfilter(params=theta,Np=5000) -> pf
      theta$loglik <- logLik(pf)
      theta
    } -> p
  p %>% arrange(Beta,mu_IR)
})-> p

p %>% 
  mutate(loglik=ifelse(loglik>max(loglik)-50,loglik,NA)) %>%
  ggplot(aes(x=Beta,y=mu_IR,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(mu[IR]))
