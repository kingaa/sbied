params <-
list(prefix = "llest/")

library(tidyverse)
library(pomp)
stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="4.2")
set.seed(1221234211)

source("https://kingaa.github.io/sbied/pfilter/model.R")

measSIR %>% pfilter(Np=1000) -> pf
logLik(pf)


if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}


bake(file="loglikest-pfilter.rds",{
  library(doRNG)
  registerDoRNG(594717807L)
  expand_grid(
    rep=1:1000,
    Np=c(1000,10000,100000)
  ) -> design
  
  foreach (p=iter(design,"row"),
           .inorder=FALSE, .combine=rbind) %dopar%
    {
      library(pomp)
      measSIR %>% pfilter(Np=p$Np) -> pf
      cbind(p,loglik=logLik(pf))
    } -> lls
  registerDoSEQ()
  lls
}) -> lls

expand_grid(
  nfilt=c(10,100,1000),
  lls
) %>%
  filter(rep<=nfilt) %>%
  ggplot(aes(x=Np,y=loglik,fill=ordered(nfilt),
             group=interaction(nfilt,Np)))+
  geom_violin(draw_quantiles=c(0.1,0.5,0.9),alpha=0.7)+
  scale_x_log10(breaks=unique(lls$Np))+
  labs(fill="nfilt")
