library(magrittr)
library(plyr)
library(readr)
library(stringi)

batchsize <- 60
replyto <- "kingaa.sismid@gmail.com"
subject <- "SISMID Module 10 Advance Instructions"
#subject <- "Trento Summer School: some advance instructions"
cmd <- 'REPLYTO=%s mutt -s "%s" -b %s -e "%s" -- %s < advance.html\n'
tweaks <- "set content_type=text/html"
  
read_csv("list.csv",col_types="cc") %>%
  mutate(
    batch=1+seq_along(email)%/%batchsize,
    firstname=stri_replace_first_regex(registrant,"(\\w+),\\s*(\\w+)","$2"),
    lastname=stri_replace_first_regex(registrant,"(\\w+),\\s*(\\w+)","$1")
    ) -> addresses

daply(
      addresses,
      ~batch,
      function (x) {
        sprintf(
          cmd,
          replyto=replyto,
          subj=subject,
          bcc=paste(x$email,collapse=","),
          tweaks=tweaks,
          to=replyto
        )
      }
) -> cmds

cat(cmds,file="spam.sh")
