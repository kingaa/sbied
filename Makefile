r-scripts: intro/intro.R stochsim/stochsim.R pfilter/pfilter.R mif/mif.R polio/polio.R measles/measles.R contacts/contacts.R ebola/ebola.R

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\",output=\"$*.R\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure

