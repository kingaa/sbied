library(pomp)

loc <- url("https://kingaa.github.io/sbied/intro/parus.csv")
dat <- read.csv(loc)
head(dat)
plot(pop~year,data=dat,type='o')

library(pomp)
parus <- pomp(dat,times="year",t0=1959)

plot(parus)

stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-c*N+e);
")
pomp(parus,rprocess=discrete_time(step.fun=stochStep,delta.t=1),
     paramnames=c("r","c","sigma"),statenames=c("N","e")) -> parus

sim <- simulate(parus,params=c(N_0=1,e_0=0,r=12,c=1,sigma=0.5),
                format="data.frame")
plot(N~year,data=sim,type='o')

rmeas <- Csnippet("pop = rpois(phi*N);")

dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")

pomp(parus,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> parus

coef(parus) <- c(N_0=1,e_0=0,r=20,c=1,sigma=0.1,phi=200)
library(ggplot2)
sims <- simulate(parus,nsim=3,format="d",include.data=TRUE)
ggplot(data=sims,mapping=aes(x=year,y=pop))+geom_line()+
  facet_wrap(~.id,ncol=1,scales="free_y")

skel <- Csnippet("DN = r*N*exp(-c*N);")

parus <- pomp(parus,skeleton=map(skel),paramnames=c("r","c"),statenames=c("N"))

traj <- trajectory(parus,params=c(N_0=1,r=12,c=1),format="data.frame")
plot(N~year,data=sim,type='o')
lines(N~year,data=traj,type='l',col='red')

plot(parus)

x <- simulate(parus)

class(x)
plot(x)

y <- as.data.frame(parus)
head(y)
head(simulate(parus,format="data.frame"))

x <- simulate(parus,nsim=10)
class(x)
sapply(x,class)
x <- simulate(parus,nsim=10,format="data.frame")
head(x)
str(x)

library(ggplot2)
x <- simulate(parus,nsim=9,format="data.frame",include.data=TRUE)
ggplot(data=x,aes(x=year,y=pop,group=.id,color=(.id=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~.id,ncol=2)

y <- trajectory(parus,format="data.frame")
plot(N~year,data=y,type="l")

coef(parus)

theta <- coef(parus)
theta[c("r","N_0")] <- c(5,3)
y <- trajectory(parus,params=theta)
plot(N~year,data=as.data.frame(y),type='l')
x <- simulate(parus,params=theta)
plot(x,var="pop")

coef(parus,c("r","N_0","sigma")) <- c(39,0.5,1)
coef(parus)
plot(simulate(parus),var="pop")
