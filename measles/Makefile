default: measles.R measles.html profile.R profile.html sigmaSE-profile1.rds sigmaSE-profile2.rds
ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

slides.pdf handout.pdf: measles.tex

.INTERMEDIATE: measles.tex

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.md\",output_format=\"html_document\")"

%.R: %.Rmd
	Rscript --vanilla -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=2)"

%.tex: %.Rnw
	Rscript -e "knitr::knit(\"$*.Rnw\",output=\"$*.tex\")"

%.pdf: %.tex
	pdflatex $*
	bibtex $*
	pdflatex $*
	pdflatex $*

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure
