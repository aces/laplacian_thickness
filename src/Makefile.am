noinst_HEADERS = \
	ParseArgv.h \
	createLaplacianGrid.h \
	laplacianGrid.h


bin_PROGRAMS = create_laplacian_grid vertex_stats laplacian_thickness

create_laplacian_grid_SOURCES = createLaplacianGrid_main.c \
	createLaplacianGrid.c 

laplacian_thickness_SOURCES = createLaplacianGrid.c laplacianGrid.cc laplacianThickness.cc

vertex_stats_SOURCES = vertexStats.cc ParseArgv.c
