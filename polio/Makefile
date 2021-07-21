default: index.html slides.pdf notes.pdf main.R \
algorithmic-parameters-exercise.html convergence-exercise.html \
demography-exercise.html initial-values-exercise.html \
starting-values-exercise.html
	make -C initial-values-exercise
	make -C starting-values-exercise

include ../rules.mk

.fresh:
	for module in initial-values-exercise starting-values-exercise; do (cd $$module && $(MAKE) fresh); done

fresh: .fresh

