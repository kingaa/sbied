params <-
list(prefix = "exercises/")

source("main.R")

measSIR %>%
  simulate(params=c(Beta=25,mu_IR=0.5,rho=0.5,k=10,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")

measSIR %>%
  simulate(params=c(Beta=40,mu_IR=0.5,rho=0.5,k=10,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")

measSIR %>%
  simulate(params=c(Beta=40,mu_IR=0.2,rho=0.5,k=10,eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")

measSIR %>%
  simulate(params=c(Beta=15,mu_IR=0.5,rho=0.5,k=10,eta=0.06,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")

seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


seir_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

measSIR %>%
  pomp(
    rprocess=euler(seir_step,delta.t=1/7),
    rinit=seir_init,
    paramnames=c("N","Beta","mu_EI","mu_IR","rho","eta"),
    statenames=c("S","E","I","R","H")
  ) -> measSEIR

measSEIR %>%
  simulate(params=c(Beta=30,mu_EI=0.8,mu_IR=1.3,rho=0.5,k=10,eta=0.06,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")

measSEIR %>% 
  simulate(params=c(Beta=40,mu_EI=0.8,mu_IR=1.3,rho=0.5,k=10,eta=0.06,N=38000),
  nsim=20,format="data.frame",include.data=TRUE) %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color="none")
