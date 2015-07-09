## ----opts,include=FALSE,cache=FALSE--------------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  pomp.cache="cache",
  encoding="UTF-8"
  )

library(ggplot2)
theme_set(theme_bw())

## ----install-pomp,eval=FALSE---------------------------------------------
## require(devtools)
## install_github("kingaa/pomp")
## # install_github("kingaa/pompExamples")

## ----prelims,cache=F-----------------------------------------------------
set.seed(594709947L)
require(ggplot2)
require(plyr)
require(reshape2)
require(pomp)
stopifnot(packageVersion("pomp")>="0.69-1")

## ----load-ricker,cache=FALSE---------------------------------------------
pompExample(ricker)

## ----plot-ricker---------------------------------------------------------
plot(ricker)

## ----sim-ricker1---------------------------------------------------------
x <- simulate(ricker)

## ------------------------------------------------------------------------
class(x)
plot(x)

## ------------------------------------------------------------------------
y <- as.data.frame(ricker)
head(y)
head(simulate(ricker,as.data.frame=TRUE))

## ------------------------------------------------------------------------
x <- simulate(ricker,nsim=10)
class(x)
sapply(x,class)
x <- simulate(ricker,nsim=10,as.data.frame=TRUE)
head(x)
str(x)

## ----fig.height=8--------------------------------------------------------
x <- simulate(ricker,nsim=9,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=x,aes(x=time,y=y,group=sim,color=(sim=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)

## ----traj-ricker---------------------------------------------------------
y <- trajectory(ricker)
dim(y)
dimnames(y)
plot(time(ricker),y["N",1,],type="l")

## ----coef-ricker---------------------------------------------------------
coef(ricker)

## ------------------------------------------------------------------------
theta <- coef(ricker)
theta[c("r","N.0")] <- c(5,3)
y <- trajectory(ricker,params=theta)
plot(time(ricker),y["N",1,],type="l")
x <- simulate(ricker,params=theta)
plot(x,var="y")

## ------------------------------------------------------------------------
coef(ricker,c("r","N.0","sigma")) <- c(39,0.5,1)
coef(ricker)
plot(simulate(ricker),var="y")

## ----bifdiag-------------------------------------------------------------
p <- parmat(coef(ricker),500)
dim(p); dimnames(p)
p["r",] <- seq(from=2,to=40,length=500)
y <- trajectory(ricker,params=p,times=200:1000)
matplot(p["r",],y["N",,],pch=".",col='black',xlab='r',ylab='N',log='x')

## ----pfilter1------------------------------------------------------------
pf <- pfilter(ricker,Np=1000)
class(pf)
plot(pf)
logLik(pf)

## ----pfilter2------------------------------------------------------------
pf <- pfilter(pf)
logLik(pf)

## ----pfilter3------------------------------------------------------------
pf <- pfilter(pf,Np=100)
logLik(pf)

## ----parus-data----------------------------------------------------------
loc <- url("http://kinglab.eeb.lsa.umich.edu/SBIED/data/parus.csv")
dat <- read.csv(loc)
head(dat)
plot(pop~year,data=dat,type='o')

## ----parus-pomp1---------------------------------------------------------
require(pomp)
parus <- pomp(dat,times="year",t0=1959)

## ----parus-plot1---------------------------------------------------------
plot(parus)

## ----parus-skel-defn-----------------------------------------------------
skel <- Csnippet("DN = r*N*exp(-N);")

## ----parus-add-skel------------------------------------------------------
parus <- pomp(parus,skeleton=skel,skeleton.type='map',
              paramnames=c("r"),statenames=c("N"))

## ----parus-first-traj,results='markup'-----------------------------------
traj <- trajectory(parus,params=c(N.0=1,r=12),as.data.frame=TRUE)
ggplot(data=traj,aes(x=time,y=N))+geom_line()

## ----parus-sim-defn------------------------------------------------------
stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-N+e);
")
pomp(parus,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
     paramnames=c("r","sigma"),statenames=c("N","e")) -> parus

## ----ricker-first-sim,results='markup'-----------------------------------
sim <- simulate(parus,params=c(N.0=1,e.0=0,r=12,sigma=0.5),
                as.data.frame=TRUE,states=TRUE)
plot(N~time,data=sim,type='o')
lines(N~time,data=traj,type='l',col='red')

## ----parus-rmeas-defn----------------------------------------------------
rmeas <- Csnippet("pop = rpois(phi*N);")

## ----parus-dmeas-defn----------------------------------------------------
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")

## ----parus-add-meas------------------------------------------------------
pomp(parus,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> parus

## ----ricker-add-params---------------------------------------------------
coef(parus) <- c(N.0=1,e.0=0,r=20,sigma=0.1,phi=200)

## ----ricker-second-sim,results='markup'----------------------------------
sims <- simulate(parus,nsim=3,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=sims,mapping=aes(x=time,y=pop))+geom_line()+
  facet_wrap(~sim)

