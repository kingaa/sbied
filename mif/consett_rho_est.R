library(tidyverse)
library(lubridate)

download.file(
  url="https://kingaa.github.io/pomp/vignettes/twentycities.rda",
  destfile=file.path(tempdir(),"twentycities.rda")
)
load(file=file.path(tempdir(),"twentycities.rda"))

measles |>
  filter(town=="Consett") -> x
demog |>
  filter(town=="Consett") -> y

x |>
  mutate(year=year(date)) |>
  group_by(year) |>
  summarize(cases=sum(cases)) |>
  left_join(y,by="year") |>
  mutate(
    cumcase=cumsum(cases),
    cumbirth=cumsum(births)
  ) -> z

z |>
  ggplot(aes(x=cumbirth,y=cumcase))+
  geom_point()+
  expand_limits(x=0,y=0)+
  geom_smooth(method="lm",formula=y~x,color="blue")+
  geom_smooth(method="lm",formula=y~x-1,color="green")

fit <- lm(cumcase~cumbirth,data=z)
summary(fit)
coef(fit)
