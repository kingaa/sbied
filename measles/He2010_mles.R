require(aakmisc)
read.csv("~/projects/measles/ms/all_estimates.csv") %>%
  ddply(~town,subset,loglik==max(loglik),
        select=-c(X,nfail.max,nfail.min,log.beta1)) %>%
  unique() %>%
  mutate(sigmaSE=1/sqrt(eta)) %>%
  rename(c(S.0="S_0",E.0="E_0",I.0="I_0",R.0="R_0",tau="psi")) %>%
  subset(select=-eta) -> mles

saveRDS(mles,file="He2010_mles.rds")

towns <- unique(mles$town)

read.csv("~/share/data/childhood/EW_Birth_Pop_1939-1964.csv") %>%
  subset(town %in% towns, select=c(town,year,pop,births)) %>%
  arrange(town,year) -> demog
read.csv("~/share/data/childhood/EW_Lat_Long.csv") %>%
  subset(town %in% towns) %>%
  rename(c(Long="long",Lat="lat")) %>%
  arrange(town) -> coord
read.csv("~/share/data/childhood/EW_Measles_1944_1964.csv",
         colClasses=c(date="Date")) %>%
  subset(town %in% towns) %>%
  subset(select=c(town,date,cases)) %>%
  arrange(town,date) %>%
  mutate(date=date+5) -> measles

save(demog,coord,measles,file="twentycities.rda",compress='xz')


