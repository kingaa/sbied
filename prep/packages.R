## check to see that the version of R is sufficiently recent
minRversion <- "3.5.1"
rv <- getRversion()
if (rv < minRversion)
  stop("R version >= ",minRversion," is required",call.=FALSE)

## get list of packages to install
pkglist <- scan(
  what=character(0),
  text="
bbmle
coda
colorspace
cowplot
deSolve
dplyr
foreach
doParallel
doRNG
ggplot2
gridExtra
gtable
knitr
lhs
magrittr
mvtnorm
nloptr
plyr
RColorBrewer
reshape2
sos
subplex
tidyr
tidyverse
pomp
"
)

## some packages may be already installed
pkglist <- setdiff(pkglist,rownames(installed.packages()))

## do installation
op <- options(warn=2)
if (length(pkglist)>0) {
  cat("trying to install packages in user directory...\n")
  lib <- Sys.getenv("R_LIBS_USER")
  dir.create(lib,recursive=TRUE,showWarnings=FALSE)
  res <- try(install.packages(pkglist,lib=lib))
  if (inherits(res,"try-error")) {
    stop("cannot install to ",lib,call.=FALSE)
  } else {
    cat("first set of packages installed successfully to user directory\n\t(",lib,")!\n")
  }
} else {
  cat("first set of packages already installed.\n")
}
options(op)
