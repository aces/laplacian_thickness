/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
#ifndef __SURFACE_FUNCS__
#define __SURFACE_FUNCS__

/***
 * functions to deal with vectors of surfaces
 * $Id: surface_functions.h,v 1.2 2009/07/27 21:48:03 claude Exp $
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
