## check to see that the version of R is at least 3.1.2
print(R.version.string)
stopifnot(getRversion()>="3.2.0")

## get list of packages to install
pkglist <- scan(
  what=character(0),
  text="
bbmle
coda
colorspace
deSolve
foreach
ggplot2
gridExtra
gtable
knitr
lhs
magrittr
maptools
mvtnorm
nloptr
plyr
RColorBrewer
reshape2
scales
sos
sp
stringr
subplex
xtable
"
)

## install the packages (you will be prompted for a mirror)
if (length(pkglist)>0) update.packages(pkglist)
install.packages(c("pomp","pompExamples"),repos="http://kinglab.eeb.lsa.umich.edu/R")

cat("all packages installed successfully!\n")
