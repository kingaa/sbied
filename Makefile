install-scripts:
	rsync -avz --chmod=a+rX,go-w scripts website

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\")"
