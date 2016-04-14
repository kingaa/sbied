## ----cores,include=FALSE,purl=TRUE,cache=FALSE---------------------------
## change this line to reflect the architecture of your machine
options(cores=30)

## ----data----------------------------------------------------------------

## ----load-package--------------------------------------------------------
library(pomp)
packageVersion("pomp")

## ----statenames----------------------------------------------------------

## ----covariates----------------------------------------------------------

## ----rp_names------------------------------------------------------------

## ----ivp_names-----------------------------------------------------------

## ----fixed_params--------------------------------------------------------

## ----param_guess---------------------------------------------------------

## ----rprocess------------------------------------------------------------

## ----measure-------------------------------------------------------------

## ----initializer---------------------------------------------------------

## ----trans---------------------------------------------------------------

## ----pomp----------------------------------------------------------------

## ----first_sim,include=FALSE---------------------------------------------
library(reshape2)
library(ggplot2)
nsim <- 9
x <- simulate(polio,nsim=nsim,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=x,mapping=aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  geom_line()+
  scale_color_manual(values=c(`TRUE`="blue",`FALSE`="red"))+
  guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)+
  scale_y_sqrt()+
  theme_bw()+theme(strip.text=element_blank()) -> pl

## ----first_pf------------------------------------------------------------
pf <- pfilter(polio,Np=1000)
logLik(pf)

## ----parallel-setup,cache=FALSE------------------------------------------
library(foreach)
library(doParallel)
registerDoParallel()

## ----pf1-----------------------------------------------------------------
set.seed(493536993,kind="L'Ecuyer")
t1 <- system.time(
  pf1 <- foreach(i=1:10,.packages='pomp',
                 .options.multicore=list(set.seed=TRUE)
  ) %dopar% {
    pfilter(polio,Np=5000)
  }
)
(L1 <- logmeanexp(sapply(pf1,logLik),se=TRUE))

## ----pf1-plot,echo=FALSE-------------------------------------------------
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
pf1 %>% 
  setNames(seq_along(pf1)) %>%
  ldply(as.data.frame,.id='rep') %>% 
  subset(select=c(time,rep,ess,cond.loglik)) %>% 
  melt(id=c('time','rep')) %>%
  ggplot(aes(x=time,y=value,group=variable))+
  geom_line()+
  facet_wrap(~variable,ncol=1,scales='free_y')+
  guides(color=FALSE)+
  theme_bw()

## ----local_search--------------------------------------------------------
stew(file="local_search.rda",{
  t1 <- system.time({
    m1 <- foreach(i=1:90,
                  .packages='pomp', .combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      mf <- mif2(polio,
                 Np=1000,
                 Nmif=50,
                 cooling.type="geometric",
                 cooling.fraction.50=0.5,
                 transform=TRUE,
                 rw.sd=rw.sd(
                   b1=0.02, b2=0.02, b3=0.02, b4=0.02, b5=0.02, b6=0.02,
                   psi=0.02, rho=0.02, tau=0.02, sigma_dem=0.02, sigma_env=0.02,
                   IO_0=ivp(0.2), SO_0=ivp(0.2)
                 )
      )
      ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  }
  )
},seed=318817883,kind="L'Ecuyer")

## ----pairs_local,fig.width=6,fig.height=6--------------------------------
pairs(~loglik+psi+rho+tau+sigma_dem+sigma_env,data=subset(m1,loglik>max(loglik)-20))

## ----param_file1---------------------------------------------------------
write.csv(m1,file="polio_params.csv",row.names=FALSE,na="")

## ----nbinom--------------------------------------------------------------
nb_lik <- function(theta) {
  -sum(dnbinom(as.vector(obs(polio)),size=exp(theta[1]),prob=exp(theta[2]),log=TRUE))
} 
nb_mle <- optim(c(0,-5),nb_lik)
-nb_mle$value

## ----arma----------------------------------------------------------------
log_y <- log(as.vector(obs(polio))+1)
arma_fit <- arima(log_y,order=c(2,0,2),seasonal=list(order=c(1,0,1),period=12))
arma_fit$loglik-sum(log_y)

## ----box-----------------------------------------------------------------
polio_box <- rbind(
  b1=c(-2,8),
  b2=c(-2,8),
  b3=c(-2,8),
  b4=c(-2,8),
  b5=c(-2,8),
  b6=c(-2,8),
  psi=c(0,0.1),
  rho=c(0,0.1),
  tau=c(0,0.1),
  sigma_dem=c(0,0.5),
  sigma_env=c(0,1),
  SO_0=c(0,1),
  IO_0=c(0,0.01)
)

## ----global_search-------------------------------------------------------
stew(file="global_search.rda",{
  t2 <- system.time({
    m2 <- foreach(i=1:400,.packages='pomp',.combine=rbind,
                  .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      guess <- apply(polio_box,1,function(x)runif(1,x[1],x[2]))
      mf <- mif2(polio,
                 start=c(guess,fixed_params),
                 Np=2000,
                 Nmif=300,
                 cooling.type="geometric",
                 cooling.fraction.50=0.5,
                 transform=TRUE,
                 rw.sd=rw.sd(
                   b1=0.02, b2=0.02, b3=0.02, b4=0.02, b5=0.02, b6=0.02,
                   psi=0.02, rho=0.02, tau=0.02, sigma_dem=0.02, sigma_env=0.02,
                   IO_0=ivp(0.2), SO_0=ivp(0.2)
                 ))
      ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
      data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
    }
  })
},seed=290860873,kind="L'Ecuyer")
library(plyr)
params <- arrange(rbind(m1,m2[names(m1)]),-loglik)
write.csv(params,file="polio_params.csv",row.names=FALSE,na="")

## ----pairs_global,fig.width=6,fig.height=6-------------------------------
pairs(~loglik+psi+rho+tau+sigma_dem+sigma_env,data=subset(m2,loglik>max(loglik)-20))

## ----plot_simulated,echo=F-----------------------------------------------
library(magrittr)
library(reshape2)
library(ggplot2)
nsim <- 9
params %>%
  subset(loglik==max(loglik)) %>%
  unlist() -> coef(polio)
polio %>%
  simulate(nsim=nsim,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  geom_line()+
  scale_color_manual(values=c(`TRUE`="blue",`FALSE`="red"))+
  guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)+
  scale_y_sqrt()+
  theme_bw()+theme(strip.text=element_blank())

## ----param_file,fig.width=6,fig.height=6---------------------------------
params <- read.csv("polio_params.csv")
pairs(~loglik+psi+rho+tau+sigma_dem+sigma_env,data=subset(params,loglik>max(loglik)-20))

## ----global_rho,echo=F---------------------------------------------------
plot(loglik~rho,data=subset(params,loglik>max(loglik)-10),log="x")

## ----load_profile_rho,include=F------------------------------------------
readRDS("profile_rho.rds") -> m3

## ----save_profile_rho----------------------------------------------------
params <- arrange(rbind(params,m3[names(params)]),-loglik)
write.csv(params,file="polio_params.csv",row.names=FALSE,na="")

## ----profile_rho_plot1,echo=F,fig.width=6,fig.height=6-------------------
params %>%
  mutate(rho.bin=cut(rho,breaks=seq(0.01,0.025,by=0.0005),include=T)) %>%
  subset(rho>0.01 & rho<0.025) %>%
  ddply(~rho.bin,subset,rank(-loglik)<=3) -> pp

# pp %>%
#   ggplot(aes(x=rho,y=loglik))+
#   geom_point()+
#   # geom_smooth(method="loess")+
#   # lims(y=max(m3$loglik)+c(-10,2))+
#   labs(x=expression(rho))+
#   theme_bw()

pairs(~loglik+psi+rho+tau+sigma_dem+sigma_env,data=subset(pp,loglik>max(loglik)-20))

## ----persistence---------------------------------------------------------
library(plyr)
library(magrittr)
library(ggplot2)
library(grid)

params %>%
  subset(loglik==max(loglik)) %>%
  unlist() -> coef(polio)

bake(file="sims.rds",seed=398906785,
     simulate(polio,nsim=2000,as.data.frame=TRUE,include.data=TRUE)
     ) -> sims
ddply(sims,~sim,summarize,zeros=sum(cases==0)) -> num_zeros

num_zeros %>%
  subset(sim != "data") %>%
  ggplot(mapping=aes(x=zeros))+
  geom_density()+
  geom_vline(data=subset(num_zeros,sim=="data"),aes(xintercept=zeros))+
  labs(x="# zero-case months")+
  theme_bw() -> pl1

num_zeros %>%
  subset(sim=="data") %>%
  extract2("zeros") -> datz

num_zeros %>%
  ddply(.(data=sim=="data"),summarize,
        mean=mean(zeros)) -> mean_zeros

sims %>%
  subset(sim != "data") %>%
  ddply(~sim,summarize,
        fadeout1=sum(IB+IO<0.5),
        fadeout80=sum(IB+IO<80)
  ) -> fadeouts

fadeouts %>%
  ggplot(mapping=aes(x=fadeout1))+
  geom_histogram(binwidth=1,fill=NA,color='black',aes(y=..density..))+
  labs(x="# fadeouts")+
  theme_bw()+theme(legend.position="top") -> pl2

print(pl1)
print(pl2,vp=viewport(x=0.8,y=0.75,height=0.4,width=0.2))

sims %>%
  subset(sim!="data") %>%
  summarize(imports=coef(polio,"psi")*mean(SO+SB1+SB2+SB3+SB4+SB5+SB6)/12
            ) -> imports

## ----check_class_m3,eval=F-----------------------------------------------
## class(m3)

## ----check_class_m3_1,eval=F---------------------------------------------
## class(m3[[1]])

## ----mif_diagnostics,eval=F----------------------------------------------
## plot(m3[r3$logLik>max(r3$logLik)-10])

## ----likelihood_convergence,eval=F---------------------------------------
## loglik_convergence <- do.call(cbind,conv.rec(m3[r3$logLik>max(r3$logLik)-10],"loglik"))
## matplot(loglik_convergence,type="l",lty=1,ylim=max(loglik_convergence,na.rm=T)+c(-10,0))

