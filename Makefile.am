INCLUDES = $(MNI_INCLUDES)
LDADD = $(MNI_LIBS) -lm

bin_PROGRAMS = create_cortical_mantle create_laplacian_volume

create_cortical_mantle_SOURCES = createMantle.cc corticalMantle.cc \
	ParseArgv.c

create_laplacian_volume_SOURCES = solve_laplace.cc create_laplace_volume.cc



