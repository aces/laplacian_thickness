bin_PROGRAMS = create_laplacian_grid laplacian_thickness_volume two_objects_evaluate vertex_stat_clean
bin_SCRIPTS = laplacian_initialisation

vertex_stat_clean_SOURCES = surfaceStats.cc ParseArgv.c

create_laplacian_grid_SOURCES = laplace_initialisation.c

laplacian_thickness_volume_SOURCES = solveLaplace.cc laplacianGrid.cc ParseArgv.c

two_objects_evaluate_SOURCES = two_objects_evaluate.c

SUFFIXES = .pl
.pl:
	@rm -f $@
	sed < $< > $@-t \
		-e 's,xPERLx,$(PERL),' -e 's,xINCDIRx,$(perllibdir),g' &&\
	mv $@-t $@

CLEANFILES = $(bin_SCRIPTS)
