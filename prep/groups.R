library(tidyverse)

read_csv("list.csv") %>%
  mutate(
    rank=rank(runif(n=n()))-1,
    ngroup=ceiling(n()/5),
    gp=(rank%%ngroup)+1,
    group=sprintf("room%02d",gp)
  ) %>%
  arrange(group,registrant) %>%
  select(
    `Pre-assign Room Name`=group,
    `Email Address`=email,
    Name=registrant
  ) -> dat

dat %>%
  select(-Name) %>%
  write_csv("zoom_rooms.csv")

dat %>%
  write_csv("groups.csv")
