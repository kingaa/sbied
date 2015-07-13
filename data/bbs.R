require(pomp)

ode <- Csnippet("                      
  int nrate = 2;
  double rate[nrate];		// transition rates
  double term[nrate];		// terms in the ODE

  // Transition rates
  rate[0] = beta*I/pop;         // force of infection
  rate[1] = gamma;		// recovery

  // compute the several terms
  term[0] = rate[1]*S;
  term[1] = rate[3]*I;

  // balance the equations
  DS = -term[0];
  DI = term[0]-term[1];
  DR = term[1];
  Dcases = term[1];		// accumulate the new I->R transitions #
")
                      
stochstep <- Csnippet("
  int nrate = 2;
  double rate[nrate];		// transition rates
  double term[nrate];		// transition numbers
  double dW;

  // Extrademographic stochasticity
  dW = rgammawn(beta_sd,dt);

  // Transition rates
  rate[0] = beta*I/pop*dW/dt;   // force of infection
  rate[1] = gamma;		// recovery

  // compute the transition numbers
  reulermultinom(1,S,&rate[0],dt,&term[0]);
  reulermultinom(1,I,&rate[1],dt,&term[1]);

  // balance the equations
  S += -term[0];
  I += term[0]-term[1];
  R += term[1];
  cases += term[1]; // accumulate the new I->R transitions #
")

dmeas <- Csnippet("
  double tol = 1e-6;
  lik = dnbinom_mu(reports,1.0/sigma/sigma,tol+rho*cases,give_log);
")

rmeas <- Csnippet("
  double tol = 1e-6;
  reports = rnbinom_mu(1.0/sigma/sigma,tol+rho*cases);
")

paruntrans <- Csnippet("
  Tgamma = log(gamma);
  Tbeta = log(beta);
  Tbeta_sd = log(beta_sd);
  Trho = logit(rho);
  Tsigma = log(sigma);
  to_log_barycentric(&TS_0,&S_0,3);
")

partrans <- Csnippet("
  Tgamma = exp(gamma);
  Tbeta = exp(beta);
  Tbeta_sd = exp(beta_sd);
  Trho = expit(rho);
  Tsigma = exp(sigma);
  from_log_barycentric(&TS_0,&S_0,3);
")

flu <- read.csv2(text="
day;reports
1;3
2;8
3;28
4;76
5;222
6;293
7;257
8;237
9;192
10;126
11;70
12;28
13;12
14;5
")

pomp(
     data=flu,
     times="day",
     t0=0,
     rprocess=euler.sim(step.fun=stochstep,delta.t=1/12),
     skeleton=ode,
     skeleton.type="vectorfield",
     rmeasure=rmeas,
     dmeasure=dmeas,
     obsnames = c("reports"),
     statenames=c("S","I","R","cases"),
     paramnames=c(
       "gamma","beta","beta.sd","pop","rho","sigma",
       "S.0","I.0","R.0"
       ),
     zeronames=c("cases"),
     parameter.inv.transform=paruntrans,
     parameter.transform=partrans,
     initializer=function(params, t0, ...) {
       fracs <-  params[c("S.0","I.0","R.0")]
       x0 <- setNames(numeric(4),c("S","I","R","cases"))
       x0[c("S","I","R")] <- round(params['pop']*fracs/sum(fracs))
       x0["cases"] <- 0
       x0
     }
     ) -> bbs

p <- c(
       gamma=1/3,beta=1.4,
       beta.sd=0,
       pop=1400,
       rho=0.9,sigma=0.01,
       S.0=0.999,I.0=0.001,R.0=0
       )

plot(cases~time,
     data=simulate(bbs,params=p,nsim=1,as.data.frame=TRUE),
     type='l')

c("bbs")
