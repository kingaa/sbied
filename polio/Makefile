pomp2E3B53FC7CC8.so:

%.so: %.c
	PKG_CFLAGS="-I `R RHOME`/library/pomp/include/" R CMD SHLIB $*.c

clean:
	$(RM) *.o *.so
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure

