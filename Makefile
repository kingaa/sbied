MODULES = prep intro stochsim pfilter mif ebola contacts measles polio misc

default: index.html syllabus.html acknowledge.html welcome.html modules README.md 

modules:
	for module in $(MODULES); do ($(MAKE) -C $$module | tee $$module.out); done

include rules.mk

welcome.html: welcome.md
	$(REXE) -e "rmarkdown::render(\"$^\",output_format=\"revealjs::revealjs_presentation\")"

clean:
	for module in $(MODULES); do $(MAKE) -C $$module clean; done

fresh:
	for module in $(MODULES); do $(MAKE) -C $$module fresh; done
