# REXE = Rscript --vanilla
# For some reason, --vanilla fails on my Mac 
REXE = Rscript --no-save --no-restore --no-init-file
RBATCH = R CMD BATCH --no-save --no-restore
PANDOC = pandoc -s -t html5+smart --mathjax
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

slides.pdf handout.pdf notes.pdf: main.tex
quiz.pdf quiz_soln.pdf: quiz_main.tex

.INTERMEDIATE: main.tex slides.tex handout.tex notes.tex \
quiz.tex quiz_soln.tex quiz_main.tex	

%.html: %.Rmd
	$(REXE) -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.md: %.Rmd
	$(REXE) -e "rmarkdown::render(\"$*.Rmd\",output_format=\"md_document\")"

%.html: %.md
	$(REXE) -e "rmarkdown::render(\"$*.md\",output_format=\"html_document\")"

%.Rout: %.R
	$(RBATCH) $*.R $@

%.Rout.save: %.R
	$(RBATCH) $*.R $@

%.R: %.Rnw
	$(REXE) -e "knitr::purl(\"$*.Rnw\",output=\"$*.R\",documentation=0)"

%.R: %.Rmd
	$(REXE) -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=0)"

%.tex: %.Rnw
	$(REXE) -e "knitr::knit(\"$*.Rnw\",output=\"$*.tex\")"

%.pdf: export BSTINPUTS=$(ROOT_DIR)

%.pdf: %.tex
	pdflatex $*
	bibtex $* || /bin/true
	pdflatex $*
	pdflatex $*

clean:
	$(RM) *.bak *.tmp
	$(RM) *.o *.so
	$(RM) *.log *.aux *.out *.blg *.toc *.nav *.snm *.vrb *.brf *.cut
	$(RM) Rplots.*

fresh: clean
	$(RM) *.bbl
	$(RM) -r tmp

