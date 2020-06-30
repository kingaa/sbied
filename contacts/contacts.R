# DEBUG <- TRUE
DEBUG <- FALSE
library(doParallel)
library(doRNG)
cores <- detectCores()
registerDoParallel(cores)
registerDoRNG(2050320976)

## contact_data <- read.table(file="contacts.csv",header=TRUE)
## matplot(t(contact_data[1:15,1:4]),
##         ylab="total sexual contacts",xlab="6-month intervals",
##         type="l",xaxp=c(1,4,3))

par(mai=c(0.8,0.8,0.1,0.1))
contact_data <- read.table(file="contacts.csv",header=TRUE)
matplot(t(contact_data[1:15,1:4]),
        ylab="total sexual contacts",xlab="6-month intervals", 
        type="l",xaxp=c(1,4,3))

knitr::include_graphics("contacts_fig4.jpg")

library(panelPomp)
contacts <- panelPompExample(pancon)

class(contacts)
slotNames(contacts)
class(unitobjects(contacts)[[1]])

coef(contacts)

timing <- !file.exists("results/pfilter1.rda")
tic <- Sys.time()

## pf1_results <- foreach(i=1:20) %dopar% {
##   library(pomp)
##   library(panelPomp)
##   pf <- pfilter(contacts,Np= if(DEBUG) 10 else 2000)
##   list(logLik=logLik(pf),
##     unitLogLik=sapply(unitobjects(pf),logLik))
## }

stew("results/pfilter1.rda",{
pf1_results <- foreach(i=1:20) %dopar% {
  library(pomp)
  library(panelPomp)
  pf <- pfilter(contacts,Np= if(DEBUG) 10 else 2000)
  list(logLik=logLik(pf),
    unitLogLik=sapply(unitobjects(pf),logLik))
}
}) 

t1 <- as.numeric(difftime(Sys.time(),tic,units="mins"))
if(timing) save(t1,file="results/pfilter1-timing.rda") else load("results/pfilter1-timing.rda")

loglik1 <- sapply(pf1_results,function(x) x$logLik)
logmeanexp(loglik1,se=T)

pf1_loglik_matrix <- sapply(pf1_results,function(x) x$unitLogLik)
panel_logmeanexp(pf1_loglik_matrix,MARGIN=1,se=T)

timing <- !file.exists("results/mif1.rda")
tic <- Sys.time()

## mif_results <- foreach(i=1:10) %dopar% {
##   library(pomp); library(panelPomp)
##   mf <- mif2(contacts,
##     Nmif = if(DEBUG) 2 else 50,
##     Np = if(DEBUG) 5 else 1000,
##     cooling.fraction.50=0.1,
##     cooling.type="geometric",
##     transform=TRUE,
##     rw.sd=rw.sd(mu_X=0.02, sigma_X=0.02,
##       mu_D = 0.02, sigma_D=0.02,
##       mu_R=0.02, sigma_R =0.02, alpha=0.02
##     )
##   )
##   list(logLik=logLik(mf),params=coef(mf))
## }

stew("results/mif1.rda",{
mif_results <- foreach(i=1:10) %dopar% {
  library(pomp); library(panelPomp)
  mf <- mif2(contacts,
    Nmif = if(DEBUG) 2 else 50,
    Np = if(DEBUG) 5 else 1000,
    cooling.fraction.50=0.1,
    cooling.type="geometric",
    transform=TRUE,
    rw.sd=rw.sd(mu_X=0.02, sigma_X=0.02,
      mu_D = 0.02, sigma_D=0.02,
      mu_R=0.02, sigma_R =0.02, alpha=0.02
    )
  )
  list(logLik=logLik(mf),params=coef(mf))
}
})

t2 <- difftime(Sys.time(),tic,units="mins")
if(timing) save(t2,file="results/mif1-timing.rda") else load("results/mif1-timing.rda")

timing <- !file.exists("results/mif1-lik-eval.rda")
tic <- Sys.time()

## mif_logLik <-  sapply(mif_results,function(x)x$logLik)
## mif_mle <- mif_results[[which.max(mif_logLik)]]$params
## pf3_loglik_matrix <- foreach(i=1:10,.combine=rbind) %dopar% {
##   library(pomp)
##   library(panelPomp)
##   unitlogLik(pfilter(contacts,shared=mif_mle,Np=if(DEBUG) 50 else 2000))
## }

stew("results/mif1-lik-eval.rda",{
mif_logLik <-  sapply(mif_results,function(x)x$logLik)
mif_mle <- mif_results[[which.max(mif_logLik)]]$params
pf3_loglik_matrix <- foreach(i=1:10,.combine=rbind) %dopar% {
  library(pomp)
  library(panelPomp)
  unitlogLik(pfilter(contacts,shared=mif_mle,Np=if(DEBUG) 50 else 2000))
}  
})

panel_logmeanexp(pf3_loglik_matrix,MARGIN=2,se=T)

t3 <- difftime(Sys.time(),tic,units="mins")
if(timing) save(t3,file="results/mif1-lik-eval-timing.rda") else load("results/mif1-lik-eval-timing.rda")
