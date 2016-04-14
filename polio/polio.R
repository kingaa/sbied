## ----data----------------------------------------------------------------
polio_data <- read.table("http://kingaa.github.io/sbied/polio/polio_wisconsin.csv")
head(polio_data)

## ----package-------------------------------------------------------------
library(pomp)
packageVersion("pomp")

## ----statenames----------------------------------------------------------
statenames <- c("SB1","SB2","SB3","SB4","SB5","SB6","IB","SO","IO")
t0 <- 1932+4/12

## ----covariates----------------------------------------------------------
bspline_basis <- periodic.bspline.basis(
  polio_data$time,nbasis=6,degree=3,period=1,names="xi%d")
covartable <- data.frame(
  time=polio_data$time,
  B=polio_data$births,
  P=predict(smooth.spline(x=1931:1954,y=polio_data$pop[12*(1:24)]),
            x=polio_data$time)$y,
  bspline_basis
)
head(covartable)

## ----rp_names------------------------------------------------------------
rp_names <- c("b1","b2","b3","b4","b5","b6","psi","rho","tau","sigma_dem","sigma_env")

## ----ivp_names-----------------------------------------------------------
ivp_names <- c("SO_0","IO_0")

## ----fixed_params--------------------------------------------------------
i <- which(abs(covartable$time-t0)<0.01)
initial_births <- as.numeric(covartable$B[i-0:5])
names(initial_births) <- c("SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0") 
fixed_params <- c(delta=1/60,initial_births)
fp_names <- c("delta","SB1_0","SB2_0","SB3_0","SB4_0","SB5_0","SB6_0")

## ----param_guess---------------------------------------------------------
params <- c(b1=3,b2=0,b3=1.5,b4=6,b5=5,b6=3,psi=0.002,rho=0.01,tau=0.001,
            sigma_dem=0.04,sigma_env=0.5,SO_0=0.12,IO_0=0.001,fixed_params)

## ----rprocess------------------------------------------------------------
rproc <- Csnippet("
  double beta = exp(dot_product(K, &xi1, &b1));
  double lambda = (beta * (IO+IB) / P + psi);
  double var_epsilon = pow(sigma_dem,2)/lambda +  sigma_env*sigma_env;
  lambda *= (var_epsilon < 1.0e-6) ? 1 : rgamma(1/var_epsilon,var_epsilon);
  double p = exp(-(delta+lambda)/12);
  double q = (1-p)*lambda/(delta+lambda);
  SB1 = B;
  SB2 = SB1*p;
  SB3 = SB2*p;
  SB4 = SB3*p;
  SB5 = SB4*p;
  SB6 = SB5*p;
  SO = (SB6+SO)*p;
  IB = (SB1+SB2+SB3+SB4+SB5+SB6)*q;
  IO = SO*q;
")

## ----measure-------------------------------------------------------------
dmeas <- Csnippet("
  double tol = 1.0e-25;
  double mean_cases = rho*IO;
  double sd_cases = sqrt(pow(tau*IO,2) + mean_cases);
  if (cases > 0.0) {
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) - pnorm(cases-0.5,mean_cases,sd_cases,1,0) + tol; 
  } else{
    lik = pnorm(cases+0.5,mean_cases,sd_cases,1,0) + tol;
  }
  if (give_log) lik = log(lik);
")

rmeas <- Csnippet("
  cases = rnorm(rho*IO, sqrt( pow(tau*IO,2) + rho*IO ) );
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")

## ----initializer---------------------------------------------------------
init <- Csnippet("
  SB1 = SB1_0;
  SB2 = SB2_0;
  SB3 = SB3_0;
  SB4 = SB4_0;
  SB5 = SB5_0;
  SB6 = SB6_0;
  IB = 0;
  IO = IO_0 * P;
  SO = SO_0 * P;
")

## ----trans---------------------------------------------------------------
toEst <- Csnippet("
 Tpsi = log(psi);
 Trho = logit(rho);
 Ttau = log(tau);
 Tsigma_dem = log(sigma_dem);
 Tsigma_env = log(sigma_env);
 TSO_0 =  logit(SO_0);
 TIO_0 = logit(IO_0);
")

fromEst <- Csnippet("
 Tpsi = exp(psi);
 Trho = expit(rho);
 Ttau = exp(tau);
 Tsigma_dem = exp(sigma_dem);
 Tsigma_env = exp(sigma_env);
 TSO_0 =  expit(SO_0);
 TIO_0 = expit(IO_0);
")

## ----pomp----------------------------------------------------------------
polio <- pomp(
  data=subset(polio_data, 
              (time > t0 + 0.01) & (time < 1953+1/12+0.01),	
              select=c("cases","time")),
  times="time",
  t0=t0,
  params=params,
  rprocess = euler.sim(step.fun = rproc, delta.t=1/12),
  rmeasure = rmeas,
  dmeasure = dmeas,
  covar=covartable,
  tcovar="time",
  statenames = statenames,
  paramnames = c(rp_names,ivp_names,fp_names),
  initializer=init,
  toEstimationScale=toEst, 
  fromEstimationScale=fromEst,
  globals="int K = 6;"
)

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

## ----pairs_local---------------------------------------------------------
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

## ----pairs_global--------------------------------------------------------
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

## ----param_file----------------------------------------------------------
params <- read.csv("polio_params.csv")
pairs(~loglik+psi+rho+tau+sigma_dem+sigma_env,data=subset(params,loglik>max(loglik)-20))

## ----global_rho----------------------------------------------------------
plot(loglik~rho,data=subset(params,loglik>max(loglik)-10),log="x")

## ----profile_rho---------------------------------------------------------
library(plyr)
library(magrittr)
library(reshape2)

bake(file="profile_rho.rda",{
  params %>% 
    subset(loglik>max(loglik)-20,select=-c(loglik,loglik.se,rho)) %>% 
    melt(id=NULL) %>% 
    daply(~variable,function(x)range(x$value)) -> box
  
  starts <- profileDesign(rho=seq(0.01,0.025,length=20),
                          lower=box[,1],upper=box[,2],
                          nprof=20)
  foreach(start=iter(starts,"row"),.packages='pomp',.combine=rbind,
          .options.multicore=list(set.seed=TRUE)
  ) %dopar% {
    mf <- mif2(polio,
               start=unlist(start),
               Np=2000,
               Nmif=100,
               cooling.type="geometric",
               cooling.fraction.50=0.1,
               transform=TRUE,
               rw.sd=rw.sd(
                 b1=0.02, b2=0.02, b3=0.02, b4=0.02, b5=0.02, b6=0.02,
                 psi=0.02, tau=0.02, sigma_dem=0.02, sigma_env=0.02,
                 IO_0=ivp(0.2), SO_0=ivp(0.2)
               ))
    mf <- mif2(mf)
    ll <- logmeanexp(replicate(10,logLik(pfilter(mf,Np=5000))),se=TRUE)
    data.frame(as.list(coef(mf)),loglik=ll[1],loglik.se=ll[2])
  }
},seed=290860873,kind="L'Ecuyer") -> m3

params <- arrange(rbind(params,m3[names(params)]),-loglik)
write.csv(params,file="polio_params.csv",row.names=FALSE,na="")

## ----profile_rho_plot1---------------------------------------------------
m3 %>% 
  ddply(.(signif(rho,3)),subset,loglik==max(loglik)) %>%
  ggplot(aes(x=rho,y=loglik))+
  geom_point()+geom_smooth()+
  # lims(y=max(m3$loglik)+c(-10,2))+
  theme_bw()

## ----persistence---------------------------------------------------------
library(plyr)
library(magrittr)
library(ggplot2)

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
  theme_bw()

num_zeros %>%
  ddply(.(data=sim=="data"),summarize,
        mean=mean(zeros)) -> mean_zeros

sims %>%
  subset(sim != "data") %>%
  ddply(~sim,summarize,
        fadeout1=sum(IB+IO<0.5),
        fadeout100=sum(IB+IO<100)
  ) -> fadeouts

fadeouts %>%
  ggplot(mapping=aes(x=fadeout1))+
  geom_histogram(binwidth=1,fill=NA,color="black",aes(y=..density..))+geom_rug()+
  theme_bw()

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

