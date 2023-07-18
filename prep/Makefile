default: index.html packages.R advance.html pre_post_tests.html

include ../rules.mk

advance.html: advance.md
	pandoc $^ -o $@

spam.sh: advance.html spam.R list.csv
	$(REXE) spam.R
