library(tidyverse)
library(stringi)

batchsize <- 30
cmd <- 'REPLYTO=%s mutt -s "%s" -b %s -e "%s" -- %s < %s\n'
replyto <- "kingaa.sismid@gmail.com"
subject <- "SISMID Module 7 Advance Instructions"
tweaks <- "set content_type=text/html"
file <- "advance.html"
## subject <- "Welcome to SISMID Module 7"
## tweaks <- ""
## file <- "welcome.md"

read_csv("list.csv",col_types="cc") |>
  mutate(
    batch=1+seq_along(email)%/%batchsize,
    firstname=stri_replace_first_regex(registrant,"(\\w+),\\s*(\\w+)","$2"),
    lastname=stri_replace_first_regex(registrant,"(\\w+),\\s*(\\w+)","$1")
  ) |>
  group_by(batch) |>
  summarize(
    command=sprintf(
      cmd,
      replyto=replyto,
      subj=subject,
      bcc=paste(email,collapse=","),
      tweaks=tweaks,
      to=replyto,
      file=file
    )
  ) |>
  pull(command) |>
  cat(file="spam.sh",sep="")
