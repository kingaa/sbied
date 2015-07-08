install-scripts:
	rsync -avz --chmod=a+rX,go-w scripts website

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure
