library(foreach)
library(doMPI)
library(pomp)

source("polio_model.R")

params <- read.csv("polio_params.csv")

cl <- startMPIcluster()
registerDoMPI(cl)

## ----profile_rho---------------------------------------------------------

library(plyr)
library(reshape2)
library(magrittr)

bake(file="profile_rho.rds",{
    params %>% 
        subset(loglik>max(loglik)-20,
               select=-c(loglik,loglik.se,rho)) %>% 
        melt(id=NULL) %>% 
        daply(~variable,function(x)range(x$value)) -> box
    
    starts <- profileDesign(rho=seq(0.01,0.025,length=30),
                            lower=box[,1],upper=box[,2],
                            nprof=10)

    foreach(start=iter(starts,"row"),
            .combine=rbind,
            .options.multicore=list(set.seed=TRUE),
            .options.mpi=list(seed=290860873,chunkSize=1)
            ) %dopar% {
                mf <- mif2(polio,
                           start=unlist(start),
                           Np=2000,
                           Nmif=300,
                           cooling.type="geometric",
                           cooling.fraction.50=0.5,
                           transform=TRUE,
                           rw.sd=rw.sd(
                               b1=0.02, b2=0.02, b3=0.02, b4=0.02, b5=0.02, b6=0.02,
                               psi=0.02, tau=0.02, sigma_dem=0.02, sigma_env=0.02,
                               IO_0=ivp(0.2), SO_0=ivp(0.2)
                           ))
                mf <- mif2(mf,Np=5000,Nmif=100,cooling.fraction.50=0.1)
                ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
                data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
            }
}) -> m3

## ----close-mpi---------------------------------------------------------

closeCluster(cl)
mpi.finalize()
