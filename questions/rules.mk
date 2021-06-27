# REXE = Rscript --vanilla
# For some reason, --vanilla fails on my Mac 
REXE = Rscript --no-save --no-restore --no-init-file

sol.pdf: sol.tex main.tex 

questions.pdf : questions.tex main.tex

.INTERMEDIATE: sol.tex main.tex questions.tex

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
	-bibtex $*
	pdflatex $*
	pdflatex $*

clean:
	$(RM) *.bak
	$(RM) *.o *.so
	$(RM) *.log *.aux *.out *.blg *.toc *.nav *.snm *.vrb *.brf
	$(RM) Rplots.*

fresh: clean
	$(RM) *.bbl
	$(RM) -r tmp

