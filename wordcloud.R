rm(list=ls())
library(tidyverse)
library(wordcloud2)
library(tm)
library(SnowballC)
library(stringr)
library(pdftools)


setwd("~/Dropbox/SISMID/Module 7/wordcloud/")

files <- dir()[grep(".pdf", dir())]

lapply(files, function(f) {
  pdf_text(f)
}) %>% unlist -> text

text %>% 
  str_replace_all("[[:punct:]]", " ") -> text

Corpus(VectorSource(text)) -> docs

# inspect(docs)

toSpace <- content_transformer(function (x, pattern) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "/")
docs <- tm_map(docs, toSpace, "@")
docs <- tm_map(docs, toSpace, "\\|")

# Convert the text to lower case
# docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove your own stop word
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)


dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
# head(d, 10)

d[d$word=="parameter",2] <- d[d$word=="parameter",2] + d[d$word=="parameters",2]
d[d$word=="rate",2] <- d[d$word=="rate",2] + d[d$word=="rates",2]
d[d$word=="model",2] <- d[d$word=="model",2] + d[d$word=="models",2]
newrows <- rbind(c("monte carlo", d[d$word=="monte",2]+d[d$word=="carlo",2]))
names(newrow) <- c("monte carlo")
d <- rbind(d, newrow)
d$freq <- as.numeric(d$freq)
d <- d %>% arrange(-freq)

## remove words manually
d <- d[!d$word %in% c("parameters", "the", "can", "using", "cases", "function", "use", "this", "results", "one", "−", 
                      "rates", "values", "error", "set", "value", "rho", "beta", "eta", "mu", "run", "units", "may",
                      "models", "mean", "ionides", "∗", "csv", "also", "lesson", "king", "small", "var", "much", "version",
                      "will", "xn−", "used", "called", "fxn", "via", "what", "following", "monte", "carlo",
                      "how", "first", "n−", "fyn", "two", "measir", "for", "library", "non", "see", "doi"),]


wordcloud2(d, minRotation = 0, maxRotation = 0, minSize = 5,
           rotateRatio = 1,color = "random-light", backgroundColor = "grey") 
