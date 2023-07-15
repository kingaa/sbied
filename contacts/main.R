library(doFuture)



par(mai=c(0.8,0.8,0.1,0.1))
contact_data <- read.table(file="contacts.csv",header=TRUE)
matplot(t(contact_data[1:15,1:4]),
  ylab="total sexual contacts",xlab="6-month intervals",
  type="l",xaxp=c(1,4,3))

library(panelPomp)
contacts <- contacts()

class(contacts)
slotNames(contacts)
class(unitobjects(contacts)[[1]])

coef(contacts)



bake("pfilter1.rds",{
  plan(multicore,workers=20)
  pf1_results <- foreach(
    i=1:20,.options.future=list(seed=TRUE)
  ) %dofuture% {
    pf <- pfilter(contacts,Np=2000)
    list(
      logLik=logLik(pf),
      unitLogLik=sapply(unitobjects(pf),logLik)
    )
  }
  attr(pf1_results,"ncpu") <- nbrOfWorkers()
  pf1_results
}) -> pf1_results
t1 <- attr(pf1_results,"system.time")[3]/60
eval_cores <- attr(pf1_results,"ncpu")

loglik1 <- sapply(pf1_results,function(x) x$logLik)
logmeanexp(loglik1,se=T)

pf1_loglik_matrix <- sapply(pf1_results,function(x) x$unitLogLik)
panel_logmeanexp(pf1_loglik_matrix,MARGIN=1,se=T)



bake("mif1.rds",{
  plan(multicore,workers=20)
  mif_results <- foreach(
    i=1:20,
    .options.future=list(seed=TRUE)
  ) %dofuture% {
    mf <- mif2(contacts,
      Nmif=50, Np=1000,
      cooling.type="geometric", # note difference with pomp
      cooling.fraction.50=0.5,
      rw.sd=rw_sd(mu_X=0.02, sigma_X=0.02,mu_D = 0.02,
        sigma_D=0.02,mu_R=0.02, alpha=0.02)
    )
    list(logLik=logLik(mf),params=coef(mf))
  }
  attr(mif_results,"ncpu") <- nbrOfWorkers()
  mif_results
}) -> mif_results
t2 <- attr(mif_results,"system.time")[3]/60
mif_cores <- attr(mif_results,"ncpu")



bake("mif1-lik-eval.rds",{
  mif_logLik <-  sapply(mif_results,function(x)x$logLik)
  mif_mle <- mif_results[[which.max(mif_logLik)]]$params
  plan(multicore,workers=10)
  pf3_loglik_matrix <- foreach(i=1:10,.combine=rbind,
    .options.future=list(seed=TRUE)
  ) %dofuture% {
    unitlogLik(pfilter(contacts,shared=mif_mle,Np=10000))
  }
  attr(pf3_loglik_matrix,"ncpu") <- nbrOfWorkers()
  pf3_loglik_matrix
}) -> pf3_loglik_matrix
t3 <- attr(pf3_loglik_matrix,"system.time")[3]/60
mif_lik_eval_cores <- attr(pf3_loglik_matrix,"ncpu")

panel_logmeanexp(pf3_loglik_matrix,MARGIN=2,se=T)
