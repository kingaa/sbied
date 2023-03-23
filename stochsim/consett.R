library(tidyverse)
library(scales)
library(lubridate)

download.file(destfile="twentycities.rda",mode='wb',
  url="http://kingaa.github.io/pomp/vignettes/twentycities.rda")
load("twentycities.rda")
file.remove("twentycities.rda")

measles |>
  group_by(town) |>
  mutate(zero=(cases==0)) |>
  summarize(n=sum(zero)) |>
  arrange(-n)

measles |>
  group_by(town) |>
  mutate(zero=(cases==0)) |>
  do({
    x <- rle(.$zero)
    tibble(len=x$lengths,val=x$values)
  }) |>
  ungroup() |>
  filter(val) |>
  select(-val) |>
  ggplot(aes(x=len))+
  geom_histogram()+
  facet_wrap(~town)

measles |>
  ggplot(aes(x=date,y=cases,group=town))+
  geom_line()+
  facet_wrap(~town)+
  scale_y_continuous(trans="log1p")

measles |>
  filter(town=="Consett",date>"1947-12-31",date<"1949-01-01") -> dat
demog |>
  filter(town=="Consett",year==1948) |>
  right_join(dat,by="town") -> dat

dat |>
  mutate(week=week(date),births=births*7/365) |>
  select(week,cases) -> dat1

dat1 |>
  write_csv(path="Measles_Consett_1948.csv")
