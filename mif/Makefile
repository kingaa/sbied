TARGETS = mif.R mif.pdf global_search.rds local_search.rds \
lik_local.rds measles_params.csv if2_settings.html 

default: $(TARGETS)
ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

.INTERMEDIATE: mif.tex

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.md\",output_format=\"html_document\")"

%.R: %.Rnw
	Rscript --vanilla -e "knitr::purl(\"$*.Rnw\",output=\"$*.R\",documentation=0)"

%.R: %.Rmd
	Rscript --vanilla -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=0)"

%.tex: %.Rnw
	Rscript --vanilla -e "knitr::knit(\"$*.Rnw\",output=\"$*.tex\")"

%.pdf: %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak *.vrb *.blg
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) *.bbl
	$(RM) -r cache figure
