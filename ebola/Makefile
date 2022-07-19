default: index.html model.html notes.pdf slides.pdf

main.tex model.html: results/forecasts.rds codes.R

results/forecasts.rds: codes.R
	Rscript codes.R
	touch $@

include ../rules.mk
