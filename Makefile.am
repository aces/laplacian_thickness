INCLUDES = $(MNI_INCLUDES)
LDADD = $(MNI_LIBS) -lm

bin_PROGRAMS = create_mantle create_laplacian_volume createMantle

createMantle_SOURCES = createMantle.cc corticalMantle.cc \
					   corticalMantle_utilities.cc

create_mantle_SOURCES = create_mantle.cc cortical_mantle.cc relax_mantle.cc \
			ParseArgv.c
create_laplacian_volume_SOURCES = solve_laplace.cc create_laplace_volume.cc



