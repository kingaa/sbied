library(magrittr)
library(plyr)
library(readr)
library(stringi)

batchsize <- 60
replyto <- "kingaa.sismid@gmail.com"
subject <- "SISMID Module 7 Advance Instructions"
cmd <- 'REPLYTO=%s mutt -s "%s" -b %s %s < msg.txt\n'

read_csv("list.csv",col_types="cc") %>%
  mutate(
    batch=1+seq_along(email)%/%batchsize,
    firstname=stri_replace_first_regex(registrant,"(\\w+),\\s+(\\w+)","$2"),
    lastname=stri_replace_first_regex(registrant,"(\\w+),\\s+\\w+","$1")
    ) -> addresses

daply(
      addresses,
      ~batch,
      function (x) {
        sprintf(cmd,replyto,subject,paste(x$email,collapse=","),replyto)
      }
      ) -> cmds

cat(cmds,file="spam.sh")
