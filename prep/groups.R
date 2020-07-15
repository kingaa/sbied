library(tidyverse)
library(stringi)

read_csv("list.csv") %>%
  mutate(
    first=stri_replace_first_regex(Registrant,"^(.+?),\\s*(.+?)$","$2") %>% stri_trans_totitle(),
    last=stri_replace_first_regex(Registrant,"^(.+?),\\s*(.+?)$","$1") %>% stri_trans_totitle(),
    slack=sprintf("@%s %s",first,last),
    rank=rank(runif(n=n()))-1,
    ngroup=ceiling(n()/5),
    gp=(rank%%ngroup)+1,
    group=sprintf("room%02d",gp)
  ) %>%
  arrange(Registrant) %>%
  select(
    Name=Registrant,
    Email="E-mail",
    Slack=slack,
    Group=group,
    Timezone
  ) -> dat

dat %>%
  select(
    `Pre-assign Room Name`=Group,
    `Email Address`=Email
  ) %>%
  write_csv("zoom_rooms.csv")

dat %>%
  write_csv("groups.csv")
