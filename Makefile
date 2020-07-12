MODULES = prep intro stochsim pfilter mif polio ebola measles contacts

default: index.html syllabus.html acknowledge.html modules

modules:
	for module in $(MODULES); do ($(MAKE) -C $$module); done

include rules.mk

.fresh:
	for module in $(MODULES); do (cd $$module && $(MAKE) fresh); done

fresh: .fresh
