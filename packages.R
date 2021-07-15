"cowplot
devtools
diagram
DiagrammeR
doParallel
doRNG
foreach
ggplot2
grid
gridExtra
iterators
knitr
lubridate
magrittr
panelPomp
pdftools
plyr
pomp
reshape2
revealjs
scales
stringi
tidyverse
" -> pkgs

readLines(textConnection(pkgs)) -> pkgs

installed.packages() |>
  rownames() -> ipkgs

pkgs |>
  setdiff(ipkgs) |>
  setdiff("") -> npkgs

install.packages(npkgs)
