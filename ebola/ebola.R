## ----prelims,include=FALSE,cache=FALSE-----------------------------------
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  pomp.cache="cache",
  encoding="UTF-8"
  )

set.seed(594709947L)
require(ggplot2)
theme_set(theme_bw())
require(plyr)
require(reshape2)
require(magrittr)
require(pomp)
stopifnot(packageVersion("pomp")>="0.69-1")

## ----get-data------------------------------------------------------------
baseurl <- "http://kinglab.eeb.lsa.umich.edu/SBIED/"
read.csv(paste0(baseurl,"data/ebola_data.csv"),stringsAsFactors=FALSE,
         colClasses=c(date="Date")) -> dat
sapply(dat,class)
head(dat)

## ----popsizes------------------------------------------------------------
## Population sizes in Guinea, Liberia, and Sierra Leone (census 2014)
populations <- c(Guinea=10628972,Liberia=4092310,SierraLeone=6190280)

## ----plot-data-----------------------------------------------------------
dat %>%
  ggplot(aes(x=date,y=cases,group=country,color=country))+
  geom_line()

## ----rproc---------------------------------------------------------------
rSim <- Csnippet('
  double lambda, beta;
  double *E = &E1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Transitions
  // From class S
  double transS = rbinom(S, 1.0 - exp(- lambda * dt)); // No of infections
  // From class E
  double transE[nstageE]; // No of transitions between classes E
  for(i = 0; i < nstageE; i++){
    transE[i] = rbinom(E[i], 1.0 - exp(-nstageE * alpha * dt));
  }
  // From class I
  double transI = rbinom(I, 1.0 - exp(-gamma * dt)); // No of transitions I->R

  // Balance the equations
  S -= transS;
  E[0] += transS - transE[0];
  for(i=1; i < nstageE; i++) {
    E[i] += transE[i-1] - transE[i];
  }
  I += transE[nstageE-1] - transI;
  R += transI;
  N_EI += transE[nstageE-1]; // No of transitions from E to I
  N_IR += transI; // No of transitions from I to R
')

## ----skel----------------------------------------------------------------
skel <- Csnippet('
  double lambda, beta;
  const double *E = &E1;
  double *DE = &DE1;
  beta = R0 * gamma; // Transmission rate
  lambda = beta * I / N; // Force of infection
  int i;

  // Balance the equations
  DS = - lambda * S;
  DE[0] = lambda * S - nstageE * alpha * E[0];
  for (i=1; i < nstageE; i++)
    DE[i] = nstageE * alpha * (E[i-1]-E[i]);
  DI = nstageE * alpha * E[nstageE-1] - gamma * I;
  DR = gamma * I;
  DN_EI = nstageE * alpha * E[nstageE-1];
  DN_IR = gamma * I;
')

## ----measmodel-----------------------------------------------------------
dObs <- Csnippet('
  double f;
  if (k > 0.0)
    f = dnbinom_mu(nearbyint(cases),1.0/k,rho*N_EI,1);
  else
    f = dpois(nearbyint(cases),rho*N_EI,1);
  lik = (give_log) ? f : exp(f);
')

rObs <- Csnippet('
  if (k > 0) {
    cases = rnbinom_mu(1.0/k,rho*N_EI);
  } else {
    cases = rpois(rho*N_EI);
  }')

## ----partrans------------------------------------------------------------
toEst <- Csnippet('
  const double *IC = &S_0;
  double *TIC = &TS_0;
  TR0 = log(R0);
  Trho = logit(rho);
  Tk = log(k);
  to_log_barycentric(TIC,IC,4);
')

fromEst <- Csnippet('
  const double *IC = &S_0;
  double *TIC = &TS_0;
  TR0 = exp(R0);
  Trho = expit(rho);
  Tk = exp(k);
  from_log_barycentric(TIC,IC,4);
')

## ----pomp-construction---------------------------------------------------
ebolaModel <- function (country=c("Guinea", "SierraLeone", "Liberia"),
                        timestep = 0.1, nstageE = 3) {

  ctry <- match.arg(country)
  pop <- unname(populations[ctry])
  nstageE <- as.integer(nstageE)

  globs <- paste0("static int nstageE = ",nstageE,";")

  dat <- subset(dat,country==ctry,select=-country)

  ## Create the pomp object
  dat %>% 
    extract(c("week","cases")) %>%
    pomp(
      times="week",
      t0=min(dat$week)-1,
      globals=globs,
      statenames=c("S","E1","I","R","N_EI","N_IR"),
      zeronames=c("N_EI","N_IR"),
      paramnames=c("N","R0","alpha","gamma","rho","k",
                   "S_0","E_0","I_0","R_0"),
      nstageE=nstageE,
      dmeasure=dObs, rmeasure=rObs,
      rprocess=discrete.time.sim(step.fun=rSim, delta.t=timestep),
      skeleton=skel, skeleton.type="vectorfield",
      toEstimationScale=toEst,
      fromEstimationScale=fromEst,
      initializer=function (params, t0, nstageE, ...) {
        all.state.names <- c("S",paste0("E",1:nstageE),"I","R","N_EI","N_IR")
        comp.names <- c("S",paste0("E",1:nstageE),"I","R")
        x0 <- setNames(numeric(length(all.state.names)),all.state.names)
        frac <- c(params["S_0"],rep(params["E_0"]/nstageE,nstageE),params["I_0"],params["R_0"])
        x0[comp.names] <- round(params["N"]*frac/sum(frac))
        x0
      }
    ) -> po
}

ebolaModel("Guinea") -> gin
ebolaModel("SierraLeone") -> sle
ebolaModel("Liberia") -> lbr

## ----load-profile--------------------------------------------------------
options(stringsAsFactors=FALSE)
profs <- read.csv(paste0(baseurl,"/ebola/ebola-profiles.csv"))

## ----profiles-plots,results='hide'---------------------------------------
require(reshape2)
require(plyr)
require(magrittr)
require(ggplot2)
theme_set(theme_bw())

profs %>% 
  melt(id=c("profile","country","loglik")) %>%
  subset(variable==profile) %>%
  ddply(~country,mutate,dll=loglik-max(loglik)) %>%
  ddply(~country+profile+value,subset,loglik==max(loglik)) %>% 
  ggplot(mapping=aes(x=value,y=dll))+
  geom_point(color='red')+
  geom_hline(yintercept=-0.5*qchisq(p=0.99,df=1))+
  facet_grid(country~profile,scales='free')+
  labs(y=expression(l))

## ----diagnostics1--------------------------------------------------------
library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
options(stringsAsFactors=FALSE)

profs %>%
  subset(country=="Guinea") %>%
  subset(loglik==max(loglik),
         select=-c(loglik,loglik.se,country,profile)) %>%
  unlist() -> coef(gin)

simulate(gin,nsim=20,as.data.frame=TRUE,include.data=TRUE) %>% 
  mutate(date=min(dat$date)+7*(time-1),
         is.data=ifelse(sim=="data","yes","no")) %>% 
  ggplot(aes(x=date,y=cases,group=sim,color=is.data,
         alpha=is.data))+
  geom_line()+
  guides(color=FALSE,alpha=FALSE)+
  scale_color_manual(values=c(no=gray(0.6),yes='red'))+
  scale_alpha_manual(values=c(no=0.5,yes=1))

## ----diagnostics-growth-rate---------------------------------------------
growth.rate <- function (y) {
  cases <- y["cases",]
  fit <- lm(log1p(cases)~seq_along(cases))
  unname(coef(fit)[2])
}
probe(gin,probes=list(r=growth.rate),nsim=500) %>% plot()

## ----diagnostics-growth-rate-and-sd--------------------------------------
growth.rate.plus <- function (y) {
  cases <- y["cases",]
  fit <- lm(log1p(cases)~seq_along(cases))
  c(r=unname(coef(fit)[2]),sd=sd(residuals(fit)))
}
probe(gin,probes=list(growth.rate.plus),
      nsim=500) %>% plot()

## ----diagnostics2,fig.height=6-------------------------------------------
log1p.detrend <- function (y) {
  cases <- y["cases",]
  y["cases",] <- as.numeric(residuals(lm(log1p(cases)~seq_along(cases))))
  y
}

probe(gin,probes=list(
  growth.rate.plus,
  probe.quantile(var="cases",prob=c(0.25,0.75)),
  probe.acf(var="cases",lags=c(1,2,3),type="correlation",
            transform=log1p.detrend)
),nsim=500) %>% plot()

## ----forecasts,eval=F----------------------------------------------------
## require(pomp)
## require(plyr)
## require(reshape2)
## require(magrittr)
## options(stringsAsFactors=FALSE)
## 
## set.seed(988077383L)
## 
## require(foreach)
## require(doMPI)
## require(iterators)
## 
## 
## horizon <- 13
## 
## foreach (country=c("SierraLeone"),.inorder=TRUE,.combine=c) %:%
##   foreach (type=c("raw","cum"),.inorder=TRUE,.combine=c) %do%
##   {
##     M1 <- ebolaModel(country=country,type=type,
##                      timestep=0.01,nstageE=3,na.rm=TRUE)
##     M2 <- ebolaModel(country=country,type="raw",
##                      timestep=0.01,nstageE=3,na.rm=TRUE)
##     time(M2) <- seq(from=1,to=max(time(M1))+horizon,by=1)
##     M3 <- ebolaModel(country=country,type="raw",
##                      timestep=0.01,nstageE=3,na.rm=TRUE)
##     time(M3) <- seq(from=max(time(M1))+1,to=max(time(M1))+horizon,by=1)
##     timezero(M3) <- max(time(M1))
##     list(M1,M2,M3)
##     } -> models
## dim(models) <- c(3,2,1)
## dimnames(models) <- list(c("fit","det.forecast","stoch.forecast"),
##                          c("raw","cum"),c("SierraLeone"))
## 
## noexport <- c("models")
## 
## ## Weighted quantile function
## wquant <- function (x, weights, probs = c(0.025,0.5,0.975)) {
##   idx <- order(x)
##   x <- x[idx]
##   weights <- weights[idx]
##   w <- cumsum(weights)/sum(weights)
##   rval <- approx(w,x,probs,rule=1)
##   rval$y
##   }
## 
## starts <- c(Guinea="2014-01-05",Liberia="2014-06-01",SierraLeone="2014-06-08")
## 
## cl <- startMPIcluster()
## registerDoMPI(cl)
## 
## bake <- function (file, expr) {
##   if (file.exists(file)) {
##     readRDS(file)
##     } else {
##       val <- eval(expr)
##       saveRDS(val,file=file)
##       val
##       }
##   }
## 
## readRDS("profiles.rds") %>%
##   ddply(~country+type+model,subset,loglik>max(loglik)-6,
##         select=-c(conv,etime,loglik.se,nfail.min,nfail.max,profile)) -> mles
## 
## mles %>% melt(id=c("country","type","model"),variable.name='parameter') %>%
##   ddply(~country+type+model+parameter,summarize,
##         min=min(value),max=max(value)) %>%
##   subset(parameter!="loglik") %>%
##   melt(measure=c("min","max")) %>%
##   acast(country~type~model~parameter~variable) -> ranges
## 
## mles %>% ddply(~country+type+model,subset,loglik==max(loglik),select=-loglik) %>%
##   mutate(k=round(k,4),rho=round(rho,4),R0=round(R0,4),E_0=3*round(E_0/3)) %>%
##   unique() %>%
##   arrange(country,type,model) -> mles
## 
## ### STOCHASTIC MODEL
## 
## bake(file="ebola-forecasts_stoch.rds",{
##   foreach (country=c("SierraLeone"),
##            .inorder=TRUE,.combine=rbind) %:%
##     foreach (type=c("raw","cum"),nsamp=c(200,200),
##              .inorder=TRUE,.combine=rbind) %do%
##     {
## 
##       params <- sobolDesign(lower=ranges[country,type,'stoch',,'min'],
##                             upper=ranges[country,type,'stoch',,'max'],
##                             nseq=nsamp)
## 
##       foreach(p=iter(params,by='row'),
##               .inorder=FALSE,
##               .combine=rbind,
##               .noexport=noexport,
##               .options.multicore=list(set.seed=TRUE),
##               .options.mpi=list(chunkSize=1,seed=1568335316L,info=TRUE)
##               ) %dopar%
##         {
##           M1 <- models["fit",type,country][[1]]
##           M2 <- models["stoch.forecast",type,country][[1]]
##           pf <- pfilter(M1,params=unlist(p),Np=2000,save.states=TRUE)
##           pf$saved.states %>% tail(1) %>% melt() %>%
##             acast(variable~rep,value.var='value') %>%
##             apply(2,function (x) {
##               setNames(c(x["S"],sum(x[c("E1","E2","E3")]),x["I"],x["R"]),
##                        c("S_0","E_0","I_0","R_0"))}) -> x
##           pp <- parmat(unlist(p),ncol(x))
##           pp[rownames(x),] <- x
##           simulate(M2,params=pp,obs=TRUE) %>%
##             melt() %>%
##             mutate(time=time(M2)[time],
##                    period=ifelse(time<=max(time(M1)),"calibration","projection"),
##                    loglik=logLik(pf))
##         } %>% subset(variable=="cases",select=-variable) %>%
##         mutate(weight=exp(loglik-mean(loglik))) %>%
##         arrange(time,rep) -> sims
## 
##       ess <- with(subset(sims,time==max(time)),weight/sum(weight))
##       ess <- 1/sum(ess^2)
##       cat("ESS stoch",country,type,"=",ess,"\n")
## 
##       sims %>% ddply(~time+period,summarize,prob=c(0.025,0.5,0.975),
##                      quantile=wquant(value,weights=weight,probs=prob)) %>%
##         mutate(prob=mapvalues(prob,from=c(0.025,0.5,0.975),
##                               to=c("lower","median","upper"))) %>%
##         dcast(period+time~prob,value.var='quantile') %>%
##         mutate(country=country,type=type)
##       }
##   }) -> fc_if
## 
## ldply(list(stoch=fc_if,det=fc_tm),.id='model') %>%
##   ddply(~country,mutate,
##         model=factor(model,levels=c("stoch","det")),
##         date=as.Date(starts[unique(as.character(country))])+7*(time-1)) %>%
##   saveRDS(file='forecasts.rds')
## 
## closeCluster(cl)
## mpi.quit()

