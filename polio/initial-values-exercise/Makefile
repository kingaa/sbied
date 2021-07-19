default: main.Rout.save cleanup

main.R: ../main.R script.patch
	cp ../main.R .
	patch < script.patch

cleanup: fresh
	$(RM) main.R.orig

include ../../rules.mk
