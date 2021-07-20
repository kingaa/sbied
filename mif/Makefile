default: index.html main.R slides.pdf notes.pdf \
if2_settings.html Q_fit_seir.html windows.html \
Q_fitall.html Q_fitall.R main_win.R
	$(RM) pomp_*.c pomp_*.so

include ../rules.mk

main_win.R: main.R win.patch
	cp main.R main_win.R
	patch -p1 main_win.R < win.patch
