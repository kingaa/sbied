library(tidyverse)
library(wordcloud2)
library(tm)
library(SnowballC)
library(stringr)
library(pdftools)

list.files(
  path=c("intro","stochsim","pfilter","mif","measles",
    "od","contacts","polio","ebola"),
  pattern=r"{.*\.pdf}",
  full.names=TRUE,
  recursive=TRUE
) -> files

files |>
  lapply(pdf_text) |>
  unlist() -> text

text |>
  Boost_tokenizer() |>
#  str_replace_all("[[:punct:]]", " ") |>
  VectorSource() |>
  Corpus() -> docs

## inspect(docs)

toSpace <- content_transformer(
  function (x, pattern) gsub(pattern, " ", x)
)

docs |>
  tm_map(toSpace, "/") |>
  tm_map(toSpace, "@") |>
  tm_map(toSpace, "\\|") |>
  ##  tm_map(content_transformer(tolower)) |>
  tm_map(removeNumbers) |>
  tm_map(removeWords, stopwords("english")) |>
  tm_map(removePunctuation) |>
  tm_map(stripWhitespace) |>
  TermDocumentMatrix() |>
  as.matrix() |>
  rowSums() |>
  sort(decreasing=TRUE) -> v

library(tidyverse)

blacklist <-
  c(
    "the", "can", "using", "cases",
    "use", "this", "one", "−",
    "values", "set", "value",
    "rho", "beta", "eta", "mu",
    "run", "may",
    "mean", "var", "loglik",
    "ionides", "based",
    "∗", "csv","np","nmif",
    "i","ii","iii","iv","v",
    "also", "lesson", "king",
    "small", "much", "version",
    "will", "xn−", "used",
    "called", "fxn", "via", "what",
    "following", 
    "how", "first", "n−",
    "fyn", "two", "measir", "for",
    "library", "librarypomp", "non", "see",
    "setrue","full",
    "doi"
  )

tibble(
  word=names(v),
  freq=v
) |>
  mutate(
    word=case_when(
      word=="parameters"~"parameter",
      word=="filter"~"filtering",
      word=="model"~"models",
      word=="monte"~"monte carlo",
      word=="carlo"~"",
      word=="rates"~"rate",
      TRUE~word
    )
  ) |>
  filter(word!="") |>
  group_by(word) |>
  summarize(freq=sum(freq)) |>
  ungroup() |>
  arrange(-freq) |>
  filter(
    ! word %in% blacklist
  ) -> dat

dat |>
  mutate(
    freq=freq/max(freq),
    freq=freq^0.8
  ) |>
  wordcloud2(
    minRotation = 0,
    maxRotation = 0,
    minSize = 1,
    rotateRatio = 1,
    color = "random-light",
    backgroundColor = "grey"
  )
