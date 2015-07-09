WEBSITE = tilia:/var/www/html/SBIED
DATAFILES = bsflu_data.txt parus.csv

install-scripts:
	rsync --delete -avz --chmod=a+rX,go-w scripts/ $(WEBSITE)/scripts

install-data:
	(cd data; rsync --delete -avz --chmod=a+rX,go-w $(DATAFILES) $(WEBSITE)/data)

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure
