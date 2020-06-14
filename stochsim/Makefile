default: stochsim.R stochsim.html bsflu_data.txt seir3_exercise.html seir3_exercise.R exer_ctmc_ode.html exercises.R exercises.html
ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript -e "rmarkdown::render(\"$*.md\",output_format=\"html_document\")"

%.R: %.Rmd
	Rscript -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=2)"

%.tex: %.Rnw
	Rscript -e "knitr::knit(\"$*.Rnw\",output=\"$*.tex\")"

%.pdf: %.tex
	pdflatex $*.tex

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure
