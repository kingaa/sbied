## TEST SCRIPT
## If this script runs without errors, pomp is usable on your system.

require(pomp)

stopifnot(packageVersion("pomp")>="0.69-1")

gomp2 <- pomp(
              data=data.frame(time=1:50,Y=NA),
              times="time",
              t0=0,
              rmeasure=Csnippet('
   Y = rlnorm(log(X),tau);
'),
              dmeasure=Csnippet('
   lik = dlnorm(Y,log(X),tau,give_log);
'),
              rprocess=discrete.time.sim(
                step.fun=Csnippet('
  double S = exp(-r*dt);
  double eps = rlnorm(0,sigma);
  X = pow(K,(1-S))*pow(X,S)*eps;
'),
                delta.t=1
                ),
              paramnames=c("sigma","tau","r","K"),
              statenames="X",
              params=c(r=0.1,K=1,sigma=0.1,tau=0.1,X.0=1)
              )

gomp2 <- simulate(gomp2)
plot(gomp2)

p <- pfilter(gomp2,Np=1000)
plot(p)

cat("'pompTest' successful!\n")
