# REXE = Rscript --vanilla
# For some reason, --vanilla fails on my Mac 
REXE = Rscript --no-save --no-restore --no-init-file

default: quiz.pdf sol.pdf 

sol.pdf: sol.tex main.tex

quiz.pdf : quiz.tex main.tex

.INTERMEDIATE: sol.tex main.tex quiz.tex comment.cut

include ../rules.mk
