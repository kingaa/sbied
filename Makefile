WEBSITE = tilia:/var/www/html/SBIED
DATAFILES = bsflu_data.txt parus.csv twentycities.rda ebola_data.csv

default: install-scripts install-data

install-scripts:
	rsync --delete -avz --chmod=a+rX,go-w scripts/ $(WEBSITE)/scripts

install-data:
	(cd data; rsync --delete -avz --chmod=a+rX,go-w $(DATAFILES) $(WEBSITE)/data)

r-scripts: intro/intro.R stochsim/stochsim.R pfilter/pfilter.R lik/lik.R mif/mif.R polio/pilio.R measles/measles.R contacts/contacts.R ebola/ebola.R

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure
