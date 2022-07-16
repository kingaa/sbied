params <-
list(prefix = "loglikest")

library(tidyverse)
library(pomp)
set.seed(1221234211)

source("https://kingaa.github.io/sbied/pfilter/model.R")

NP <- 10000
NREPS <- 10
timer <- system.time(
  pf <- replicate(
    NREPS,
    measSIR %>% pfilter(Np=NP)
  )
)
ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)





bake(file="loglikest-pfilter.rds",{
  if (file.exists("CLUSTER.R")) {
    source("CLUSTER.R")
  } else {
    library(doParallel)
    registerDoParallel()
  }
  library(doRNG)
  registerDoRNG(594717807L)
  expand_grid(
    rep=1:1000,
    NP=c(1000,10000,100000)
  ) -> design
  
  foreach (p=iter(design,"row"),
           .inorder=FALSE, .combine=rbind) %dopar%
    {
      library(pomp)
      measSIR %>% pfilter(Np=p$NP) -> pf
      cbind(p,loglik=logLik(pf))
    } -> lls
  registerDoSEQ()
  lls
}) -> lls

expand_grid(
  NREPS=c(10,100,1000),
  lls
) %>%
  filter(rep<=NREPS) %>%
  ggplot(aes(x=NP,y=loglik,fill=ordered(NREPS),
             group=interaction(NREPS,NP)))+
  geom_violin(draw_quantiles=c(0.1,0.5,0.9),alpha=0.7)+
  scale_x_log10(breaks=unique(lls$NP))+
  labs(fill="NREPS",x="NP")



sd1 <- replicate(10000,logmeanexp(rnorm(10,mean=0,sd=1),se=TRUE))
sd5 <- replicate(10000,logmeanexp(rnorm(10,mean=0,sd=5),se=TRUE))
m1 <- mean(sd1[1,])
t1 <- (sd1[1,]-m1)/sd1[2,]
m5 <- mean(sd5[1,])
t5 <- (sd5[1,]-m5)/sd5[2,]
x_range <- range(c(t1,t5))
par(mfrow=c(2,1))
par(mai=c(0.5,0.5,0.5,0.5))
hist(t1,breaks=50,xlim=x_range,main="Error in SE units, with sd=1")
abline(v=c(-2,2),col="red")
hist(t5,breaks=50,xlim=x_range,main="Error in SE units, with sd=5")
abline(v=c(-2,2),col="red")
