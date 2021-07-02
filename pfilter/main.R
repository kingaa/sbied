library(tidyverse)
library(pomp)
options(stringsAsFactors=FALSE)
stopifnot(packageVersion("pomp")>="3.4")
set.seed(1350254336)












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
  slice_design(
    center=coef(measSIR),
    Beta=rep(seq(from=5,to=30,length=40),each=3),
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
  guides(color="none")+
  labs(x="parameter value",color="")







bake(file="pfilter-grid1.rds",{
  expand.grid(
    Beta=rep(seq(from=10,to=30,length=40),each=3),
    mu_IR=rep(seq(from=0.4,to=1.5,length=40),each=3),
    rho=0.5,k=10,eta=0.06,N=38000
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
  mutate(loglik=ifelse(loglik>max(loglik)-25,loglik,NA)) %>%
  ggplot(aes(x=Beta,y=mu_IR,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(mu[IR]))
