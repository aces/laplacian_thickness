#ifndef _RELAX_MANTLE_H_
#define _RELAX_MANTLE_H_

extern "C" {
#include "volume_io.h"
#include "bicpl.h"
}

#include <iostream.h>
#include "surface_volume_defs.h"

int determine_neighbouring_mantle_points(Volume cortical_mantle,
					 int *sizes);
void relax_mantle(Volume cortical_mantle, 
		  int *sizes, 
		  int max_iterations,
		  int convergence_criteria);

#endif // _RELAX_MANTLE_H_
