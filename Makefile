default: r-scripts html-docs datafiles results

r-scripts: intro/intro.R stochsim/stochsim.R pfilter/pfilter.R mif/mif.R polio/polio.R measles/measles.R contacts/contacts.R ebola/ebola.R measles/measles-profile.R
	for i in $^; do cp $$i www/$$i; done

html-docs: index.html intro/intro.html stochsim/stochsim.html pfilter/pfilter.html mif/mif.html polio/polio.html measles/measles.html contacts/contacts.html ebola/ebola.html pfilter/monteCarlo.html measles/measles-profile.html prep/preparation.html
	for i in $^; do cp $$i www/$$i; done

datafiles: data/ebola_data.csv data/parus.csv contacts/contacts.csv contacts/contacts2.csv polio/polio_wisconsin.csv 
	for i in $^; do cp $$i www/$$i; done

results: contacts/mif1.rda contacts/pfilter1.rda contacts/pfilter2.rda mif/box_search_global.rda mif/box_search_local.rda mif/lik_local.rda mif/pf.rda polio/local_search.rda polio/global_search.rda polio/sims.rds measles/sigmaSE-profile1.rds measles/sigmaSE-profile2.rds polio/polio_params.csv ebola/ebola-profiles.csv mif/bsflu_params.csv
	for i in $^; do cp $$i www/$$i; done

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.md\",output_format=\"html_document\")"

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\",output=\"$*.R\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	find . -name "cache" -print | xargs rm -rf
	find . -name "figure" -print | xargs rm -rf
