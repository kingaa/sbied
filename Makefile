MODULES = prep misc intro stochsim pfilter mif ebola contacts measles polio

default: index.html acknowledge.html modules README.md 

modules:
	for module in $(MODULES); do ($(MAKE) -C $$module | tee $$module.out); done

serve:
	jekyll serve -b /sbied

include rules.mk

welcome.html: welcome.md
	$(REXE) -e "rmarkdown::render(\"$^\",output_format=\"revealjs::revealjs_presentation\")"

clean: .clean
	for module in $(MODULES); do $(MAKE) -C $$module clean; done

fresh: .fresh
	for module in $(MODULES); do $(MAKE) -C $$module fresh; done
