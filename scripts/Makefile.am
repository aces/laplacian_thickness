bin_SCRIPTS = laplacian_initialisation


SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' -e 's,xINCDIRx,$(perllibdir),g' &&\
	mv $@-t $@

CLEANFILES = $(bin_SCRIPTS)

