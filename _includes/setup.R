stopifnot(getRversion()>="4.0")
stopifnot(packageVersion("pomp")>="4.2")

library(knitr)
if (!exists("params")) params <- list()
opts_chunk$set(
  cache=TRUE,
  cache.path=paste("tmp",as.character(params$prefix),"cache","",sep="/"),
  comment=NA,
  echo=TRUE,
  eval=TRUE,
  dev="png",
  dev.args=list(bg="transparent"),
  dpi=300,
  error=FALSE,
  fig.align="center",
  fig.height=4,fig.width=6.83,
  fig.lp="fig:",
  fig.path=paste("tmp",as.character(params$prefix),"figure","",sep="/"),
  fig.pos="h!",
  fig.show="asis",
  highlight=TRUE,
  include=TRUE,
  message=FALSE,
  progress=TRUE,
  prompt=FALSE,
  purl=TRUE,
  results="markup",
  size="small",
  strip.white=TRUE,
  tidy=FALSE,
  warning=FALSE
  )

options(
  width=60, # number of characters in R output before wrapping
  keep.source=TRUE,
  encoding="UTF-8",
  pomp_archive_dir="results"
)

library(ggplot2)
theme_set(theme_bw())

knit_hooks$set(
  document = function (x) {
    sub("\\usepackage[]{color}","\\usepackage{xcolor}",x,fixed=TRUE)
  }
)

registerS3method(
  "knit_print",
  "data.frame",
  function (x, ...) {
    print(x,row.names=FALSE)
  }
)

myround<- function (x, digits = 1) {
  # adapted from the broman package
  # solves the bug that round() kills significant trailing zeros
  if (length(digits) > 1) {
    digits <- digits[1]
    warning("Using only digits[1]")
  }
  if (digits < 1) {
    as.character(round(x,digits))
  } else {
    tmp <- sprintf(paste("%.", digits, "f", sep = ""), x)
    zero <- paste0("0.", paste(rep("0", digits), collapse = ""))
    tmp[tmp == paste0("-", zero)] <- zero
    tmp
  }
}

mysignif <- function (x, digits = 1) {
  myround(x, digits - ceiling(log10(abs(x))))
}
