#ifndef __SURFACE_FUNCS__
#define __SURFACE_FUNCS__

/***
 * functions to deal with vectors of surfaces
 * $Id$
 ***/

using namespace std;

extern "C" {
#include <bicpl.h>
}

#include <iostream>
#include <vector>
#include <fstream>

extern void resample_polygon_vector( vector<polygons_struct> &layers,
				     unsigned short subsampling_steps );
extern void create_polygon_bintrees( vector<polygons_struct> &layers );
extern void create_intermediate_surfaces( polygons_struct *surface1,
					  polygons_struct *surface2,
					  unsigned short num_layers,
					  vector<polygons_struct> &layers );


#endif // __SURFACE_FUNCTIONS__
