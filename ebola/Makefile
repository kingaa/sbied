default: index.html model.html notes.pdf slides.pdf

main.tex model.html: codes.R forecasts.rds

forecasts.rds: codes.R
	Rscript codes.R

include ../rules.mk
