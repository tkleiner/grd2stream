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


maintainer-check:
	./configure --enable-gmt-api  && $(MAKE) clean && $(MAKE)
	./configure --disable-gmt-api && $(MAKE) clean && $(MAKE)

# # Get git commit version and date
# # found on
# # https://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automatically-print-out-its-git-version-hash
# GIT_VERSION := $(shell git --no-pager describe --tags --always --dirty)
# GIT_DATE := $(firstword $(shell git --no-pager show --date=short --format="%ad" --name-only))
#
#  # recompile version.h dependants when GIT_VERSION changes, uses temporary file version~
# .PHONY: force
# git_version~: force
# 	@echo '$(GIT_VERSION) $(GIT_DATE)' | cmp -s - $@ || echo '$(GIT_VERSION) $(GIT_DATE)' > $@
# git_version.h: git_version~
# 	@touch $@
# 	@echo Git version $(GIT_VERSION) $(GIT_DATE)
