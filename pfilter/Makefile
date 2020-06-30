TARGETS = slides.pdf handout.pdf notes.pdf pfilter.R \
bootstrap.html deterministic.html \
directSimulation.html directSimulation.R fisherSE.html \
monteCarlo.html pf-in-Nelder-Mead.html pf-in-Nelder-Mead.R \
loglikest.html loglikest.R expense.R basic_exercise_pfilter.html \
basic_exercise_pfilter.R
REXE = Rscript --vanilla

default: $(TARGETS)
ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

slides.pdf handout.pdf notes.pdf: pfilter.tex

.INTERMEDIATE: pfilter.tex

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	$(REXE) -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	$(REXE) -e "rmarkdown::render(\"$*.md\",output_format=\"html_document\")"

%.R: %.Rnw
	$(REXE) -e "knitr::purl(\"$*.Rnw\",output=\"$*.R\",documentation=0)"

%.R: %.Rmd
	$(REXE) -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=0)"

%.tex: %.Rnw
	$(REXE) -e "knitr::knit(\"$*.Rnw\",output=\"$*.tex\")"

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
