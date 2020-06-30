TARGETS= main.R notes.pdf slides.pdf handout.pdf profile.R profile.html \
sigmaSE-profile1.rds sigmaSE-profile2.rds
REXE = Rscript --vanilla

default: $(TARGETS)
ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

slides.pdf handout.pdf notes.pdf: main.tex

.INTERMEDIATE: main.tex notes.tex slides.tex handout.tex

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
