params <-
list(prefix = "seir3_exercise/")

library(tidyverse)
library(pomp)
stopifnot(packageVersion("pomp")>="3.0")
theme_set(theme_bw())
options(stringsAsFactors=FALSE)
set.seed(1221234211)



library(tidyverse)
library(pomp)

bsflu %>%
  select(day,B) %>%
  pomp(times="day",t0=-6,
    rprocess=euler(Csnippet("
      double N = 763;
      double t1 = rbinom(S,1-exp(-Beta*I/N*dt));
      double t2 = rbinom(E,1-exp(-mu_E*dt));
      double t3 = rbinom(I,1-exp(-mu_I*dt));
      double t4 = rbinom(R1,1-exp(-mu_R1*dt));
      double t5 = rbinom(R2,1-exp(-mu_R2*dt));
      S  -= t1;
      E  += t1 - t2;
      I  += t2 - t3;
      R1 += t3 - t4;
      R2 += t4 - t5;"),
      delta.t=1/5),
    rinit=Csnippet("
      S  = 762;
      E  = 0;
      I  = 1;
      R1 = 0;
      R2 = 0;"),
    dmeasure=Csnippet("
      lik = dpois(B,rho*R1+1e-6,give_log);"),
    rmeasure=Csnippet("
      B = rpois(rho*R1+1e-6);"),
    statenames=c("S","E","I","R1","R2"),
    paramnames=c("Beta","mu_E","mu_I","mu_R1","mu_R2","rho")
  ) -> flu

with(bsflu,c(mu_R1=1/(sum(B)/512),mu_R2=1/(sum(C)/512)))

coef(flu) <- c(Beta=5,mu_E=0.5,mu_I=1,mu_R1=0.33,mu_R2=0.55,rho=0.95)

flu %>%
  simulate(nsim=20,format="data.frame",include.data=TRUE) -> simdat

simdat %>%
  select(day,B,.id) %>%
  ggplot(aes(x=day,y=B,color=(.id=="data"),size=(.id=="data"),group=.id))+
  geom_line()+
  scale_color_manual(values=c(`TRUE`="red",`FALSE`="black"))+
  scale_size_manual(values=c(`TRUE`=1,`FALSE`=0.5))+
  guides(color=FALSE,size=FALSE)
