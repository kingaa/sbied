params <-
list(prefix = "basic_exercise/")

set.seed(594709947L)
library(ggplot2)
library(plyr)
library(reshape2)
library(magrittr)
library(pomp)
stopifnot(packageVersion("pomp")>="3.0")

source("https://kingaa.github.io/sbied/pfilter/model.R")

NP <- 50000
REPLICATES <- 10
timer <- system.time(
  pf <- replicate(
    REPLICATES,
    measSIR %>% pfilter(Np=NP)
  )
)
ll <- sapply(pf,logLik)
logmeanexp(ll,se=TRUE)

logmeanexp

set.seed(23)
sd5 <- replicate(10000,logmeanexp(rnorm(10,mean=0,sd=5),se=TRUE))
sd1 <- replicate(10000,logmeanexp(rnorm(10,mean=0,sd=1),se=TRUE))
m5 <- mean(sd5[1,])
t5 <- (sd5[1,]-m5)/sd5[2,]
m1 <- mean(sd1[1,])
t1 <- (sd1[1,]-m1)/sd1[2,]
x_range <- range(c(t1,t5))
par(mfrow=c(2,1))
par(mai=c(0.5,0.5,0.5,0.5))
hist(t5,breaks=50,xlim=x_range,main="Error in SE units, with sd=5")
abline(v=c(-2,2),col="red")
hist(t1,breaks=50,xlim=x_range,main="Error in SE units, with sd=1")
abline(v=c(-2,2),col="red")
