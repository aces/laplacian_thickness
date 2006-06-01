#ifndef __CREATE_LAPLACIAN_GRID__
#define __CREATE_LAPLACIAN_GRID__

#include  <volume_io/internal_volume_io.h>
#include  <bicpl.h>

extern int create_mantle ( char *input_volume_filename,
			   char *grey_surface_filename,
			   char *white_surface_filename,
			   int include_white_boundary,
			   int include_grey_boundary,
                           Volume vol );

#endif

