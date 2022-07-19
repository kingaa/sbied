MODULES = prep intro stochsim pfilter mif ebola contacts measles polio misc

default: index.html syllabus.html acknowledge.html welcome.html modules README.md 

modules:
	for module in $(MODULES); do ($(MAKE) -C $$module); done

include rules.mk

.fresh:
	for module in $(MODULES); do (cd $$module && $(MAKE) fresh); done

fresh: .fresh

welcome.html: welcome.md
	$(REXE) -e "rmarkdown::render(\"$^\",output_format=\"revealjs::revealjs_presentation\")"
