DEBUG <- FALSE
# DEBUG <- TRUE
library(doParallel)
library(doRNG)
cores <- detectCores()
registerDoParallel(cores)
registerDoRNG(2050320976)



par(mai=c(0.8,0.8,0.1,0.1))
contact_data <- read.table(file="contacts.csv",header=TRUE)
matplot(t(contact_data[1:15,1:4]),
        ylab="total sexual contacts",xlab="6-month intervals", 
        type="l",xaxp=c(1,4,3))



library(panelPomp)
contacts <- panelPompExample(pancon)

class(contacts)
slotNames(contacts)
class(unitobjects(contacts)[[1]])

coef(contacts)



stew("results/pfilter1.rda",{
  tic <- Sys.time()
  eval_cores <- cores
  pf1_results <- foreach(i=1:20) %dopar% {
    library(panelPomp)
    pf <- pfilter(contacts,Np= if(DEBUG) 10 else 2000)
    list(logLik=logLik(pf),
         unitLogLik=sapply(unitobjects(pf),logLik))
  }
  t1 <- as.numeric(difftime(Sys.time(),tic,units="mins"))
}) 

loglik1 <- sapply(pf1_results,function(x) x$logLik)
logmeanexp(loglik1,se=T)

pf1_loglik_matrix <- sapply(pf1_results,function(x) x$unitLogLik)
panel_logmeanexp(pf1_loglik_matrix,MARGIN=1,se=T)



stew("results/mif1.rda",{
  mif_cores <- cores
  tic <- Sys.time()
mif_results <- foreach(i=1:20) %dopar% {
  library(pomp); library(panelPomp)
  mf <- mif2(contacts,
    Nmif = if(DEBUG) 2 else 50,
    Np = if(DEBUG) 5 else 1000,
    cooling.type="geometric", # needed for panelPomp 0.10
    cooling.fraction.50=0.5,
    rw.sd=rw.sd(mu_X=0.02, sigma_X=0.02,mu_D = 0.02,
      sigma_D=0.02,mu_R=0.02, alpha=0.02)
  )
  list(logLik=logLik(mf),params=coef(mf))
}
  t2 <- difftime(Sys.time(),tic,units="mins")
})



stew("results/mif1-lik-eval.rda",{
tic <- Sys.time()
mif_logLik <-  sapply(mif_results,function(x)x$logLik)
mif_mle <- mif_results[[which.max(mif_logLik)]]$params
pf3_loglik_matrix <- foreach(i=1:10,.combine=rbind) %dopar% {
  library(panelPomp)
  unitlogLik(pfilter(contacts,
    shared=mif_mle,Np=if(DEBUG) 50 else 10000))
}  
t3 <- difftime(Sys.time(),tic,units="mins")
mif_lik_eval_cores <- cores
})

panel_logmeanexp(pf3_loglik_matrix,MARGIN=2,se=T)
