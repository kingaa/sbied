## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  pomp.cache="cache",
  encoding="UTF-8"
  )

## ----prelims,echo=F,cache=F----------------------------------------------
set.seed(594709947L)
require(ggplot2)
theme_set(theme_bw())
require(plyr)
require(reshape2)
require(foreach)
require(doMC)
require(pomp)
stopifnot(packageVersion("pomp")>="0.69-1")

## ----flu-data1-----------------------------------------------------------
baseurl <- "http://kinglab.eeb.lsa.umich.edu/SBIED/"
url <- paste0(baseurl,"data/bsflu_data.txt")
bsflu <- read.table(url)
head(bsflu)

## ----flu-data2-----------------------------------------------------------
ggplot(data=bsflu,aes(x=day,y=B))+geom_line()+geom_point()

## ----sir-diagram,echo=FALSE,cache=FALSE----------------------------------
require(DiagrammeR)
DiagrammeR("graph LR; S(S) --> I; I(I) --> R(R);"
           ,height=200,width=500)

## ----rproc1--------------------------------------------------------------
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
")

## ----init1---------------------------------------------------------------
sir_init <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
")

## ----rproc1-pomp---------------------------------------------------------
pomp(bsflu,time="day",t0=0,rprocess=euler.sim(sir_step,delta.t=1/6),
     initializer=sir_init,paramnames=c("N","beta","gamma"),
     statenames=c("S","I","R")) -> sir

## ----rproc2--------------------------------------------------------------
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = N-1;
  I = 1;
  R = 0;
  H = 0;
")

pomp(sir,rprocess=euler.sim(sir_step,delta.t=1/6),initializer=sir_init,
     paramnames=c("beta","gamma","N"),statenames=c("S","I","R","H")) -> sir

## ----zero1---------------------------------------------------------------
pomp(sir,zeronames="H") -> sir

## ----meas-model----------------------------------------------------------
dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
rmeas <- Csnippet("B = rbinom(H,rho);")

## ----add-meas-model------------------------------------------------------
sir <- pomp(sir,rmeasure=rmeas,dmeasure=dmeas,statenames="H",paramnames="rho")

## ------------------------------------------------------------------------
sims <- simulate(sir,params=c(beta=1.5,gamma=1,rho=0.9,N=2600),
                 nsim=20,as=TRUE,include=TRUE)

ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

## ----seir-diagram,echo=FALSE,cache=FALSE---------------------------------
require(DiagrammeR)
DiagrammeR("graph LR; S(S) --> E; E(E) --> I; I(I) --> R(R);"
           ,height=200,width=600)

## ----bsflu-plot2---------------------------------------------------------
require(reshape2)
ggplot(data=melt(bsflu,id="day"),mapping=aes(x=day,y=value,color=variable))+
  geom_line()+geom_point()

## ----sirr-diagram,echo=FALSE,cache=FALSE---------------------------------
require(DiagrammeR)
DiagrammeR("graph LR; S(S) --> I; I(I) --> R1(R1); R1 --> R2(R2);"
           ,height=200,width=600)

