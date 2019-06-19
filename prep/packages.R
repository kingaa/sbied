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
scales
sos
subplex
tidyr
tidyverse
pomp
"
)

lib <- Sys.getenv("R_LIBS_USER")

inst_pkg <- function (pkglist, lib = Sys.getenv("R_LIBS_USER")) {
  op <- options(warn=2)

  pkglist <- setdiff(pkglist,rownames(installed.packages()))
  
  if (length(pkglist)>0) {
    cat("trying to install packages in user directory...\n")
    dir.create(lib,recursive=TRUE,showWarnings=FALSE)
    res <- try(install.packages(pkglist,lib=lib))
    if (inherits(res,"try-error")) {
      stop("cannot install to ",lib,call.=FALSE)
    }
  }

  options(op)
  invisible(NULL)
}

## latest 'rngtools' requires version >= 3.6.0
if (rv < "3.6.0") {
  
  inst_pkg(c("pkgmaker","stringr","digest"),lib=lib)
  
  pkgurl <- "https://cran.r-project.org/src/contrib/Archive/rngtools/rngtools_1.3.1.tar.gz"
  install.packages(pkgurl,repos=NULL,type="source",lib=lib)
  
}

inst_pkg(pkglist,lib=lib)
cat("first set of packages installed successfully to user directory\n\t(",lib,")!\n")
