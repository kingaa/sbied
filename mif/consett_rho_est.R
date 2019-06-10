library(tidyverse)
library(lubridate)

load("~/projects/Rpkg/pomp/www/vignettes/twentycities.rda")

measles %>%
  filter(town=="Consett") -> x
demog %>%
  filter(town=="Consett") -> y

x %>%
  mutate(year=year(date)) %>%
  group_by(year) %>%
  summarize(cases=sum(cases)) %>%
  left_join(y,by="year") %>%
  mutate(
    cumcase=cumsum(cases),
    cumbirth=cumsum(births)
  ) -> z

plot(cumcase~cumbirth,data=z)
fit <- lm(cumcase~cumbirth,data=z)
summary(fit)
coef(fit)
