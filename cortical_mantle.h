#ifndef _CORTICAL_MANTLE_H_
#define _CORTICAL_MANTLE_H_

// standard headers
#include <iostream.h>
#include <map>
#include <list>

// minc headers
extern "C" {
#include "volume_io.h"
#include "bicpl.h"
}

// my headers
#include "relax_mantle.h"
#include "surface_volume_defs.h"
#include "jplPoint.h"
#include "jplPointMap.h"

int check_if_between_surfaces(Volume object_intersects,
			      int* sizes, 
			      int coord1, int coord2, int coord3, 
			      int longest_line);

int expand_from_point_in_mantle(Volume object_intersects, 
				Volume cortical_mantle,
				int* sizes,
				int coord1, int coord2, int coord3);

int find_cortical_mantle(Volume object_intersects);

#endif // _CORTICAL_MANTLE_

