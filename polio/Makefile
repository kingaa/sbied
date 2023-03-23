default: index.html slides.pdf notes.pdf main.R \
algorithmic-parameters-exercise.html convergence-exercise.html \
demography-exercise.html initial-values starting-values \
runlevel1 runlevel2
	make -C runlevel1
	make -C runlevel2

export RUNLEVEL = 3

include ../rules.mk

starting-values:
	make -C starting-values-exercise
	make starting-values-exercise.html

initial-values:
	make -C initial-values-exercise
	make initial-values-exercise.html

fresh: .fresh
	for module in initial-values-exercise starting-values-exercise runlevel1 runlevel2; do (cd $$module && $(MAKE) fresh); done
