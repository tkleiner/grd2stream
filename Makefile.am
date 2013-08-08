##
## Time-stamp: <2013-08-08 09:32:35 (tkleiner)>
##
ACLOCAL_AMFLAGS  = -I m4
AUTOMAKE_OPTIONS = foreign

SUBDIRS = m4 src example


##
## However maintainer-clean should not delete 
## anything that needs to exist in order to 
## run ‘./configure && make’. 
##
## Remove all auto-generated files
##
.PHONY: total-clean
total-clean: maintainer-clean
	@echo "make total-clean is intended for maintainers to use"
	@echo "it deletes files that require bootstrapping to rebuild."
	-rm -f aclocal.m4 config.h.in configure
	-rm -rf admin
	find -name Makefile.in | xargs rm -f
	@echo "run './bootstrap.sh' to generate files using gnu autotools"	
