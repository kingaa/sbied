TARGETS = slides.pdf handout.pdf notes.pdf main.R \
exercises.html exercises.R exer_ctmc_ode.html \
seir3_exercise.html seir3_exercise.R
REXE = Rscript

default: $(TARGETS)
ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

slides.pdf handout.pdf notes.pdf: main.tex

.INTERMEDIATE: main.tex slides.tex handout.tex notes.tex

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
	$(RM) -r tmp
