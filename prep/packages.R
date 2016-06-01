## check to see that the version of R is at least 3.2.1
print(R.version.string)
res <- try(stopifnot(getRversion()>="3.2.1"))
if (inherits(res,"try-error")) {
    cat("Please install a more recent version of R!\n")
} else {
    cat("Your version of R is sufficiently recent.\n")
}    

## get list of packages to install
pkglist <- scan(
  what=character(0),
  text="
bbmle
coda
colorspace
deSolve
foreach
doParallel
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
stringr
subplex
xtable
"
)

pkglist <- setdiff(pkglist,rownames(installed.packages()))

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
