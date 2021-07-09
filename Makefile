MODULES = prep intro stochsim pfilter mif measles polio ebola contacts quiz

default: index.html syllabus.html acknowledge.html welcome.html modules

modules:
	for module in $(MODULES); do ($(MAKE) -C $$module); done

include rules.mk

.fresh:
	for module in $(MODULES); do (cd $$module && $(MAKE) fresh); done

fresh: .fresh

welcome.html: welcome.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	$(REXE) -e "rmarkdown::render(\"$^\",output_format=\"revealjs::revealjs_presentation\")"
