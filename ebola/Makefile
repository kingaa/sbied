default: index.html model.html notes.pdf slides.pdf

main.tex model.html: forecasts.rds codes.R

forecasts.rds: codes.R
	Rscript codes.R
	touch $@

include ../rules.mk
