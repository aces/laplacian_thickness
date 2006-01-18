bin_SCRIPTS = laplacian_initialisation

EXTRA_DIST = laplacian_initialisation.pl

SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' -e 's,xINCDIRx,$(perllibdir),g' &&\
	mv $@-t $@

CLEANFILES = $(bin_SCRIPTS)

