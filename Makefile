default: r-scripts html-docs

r-scripts: intro/intro.R stochsim/stochsim.R pfilter/pfilter.R mif/mif.R polio/polio.R measles/measles.R contacts/contacts.R ebola/ebola.R measles/measles-profile.R

html-docs: index.html intro/intro.html stochsim/stochsim.html pfilter/pfilter.html mif/mif.html polio/polio.html measles/measles.html contacts/contacts.html ebola/ebola.html pfilter/monteCarlo.html measles/measles-profile.html

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\",output=\"$*.R\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	find . -name "cache" -print | xargs rm -rf
	find . -name "figure" -print | xargs rm -rf
