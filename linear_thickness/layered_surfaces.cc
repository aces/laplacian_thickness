/****
 * Code in this file deals with creating a series of layered surfaces.
 *
 * $Id$
 ****/

extern "C" {
#include <bicpl.h>
#include <ParseArgv.h>
}
#include "surface_functions.h"

#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

// argument parsing defaults
int num_subdivisions        = 2;
int num_layers              = 2;

// argument parsing table
ArgvInfo argTable[] = {
  { "-num_subdivisions", ARGV_INT, (char*)0, (char*) &num_subdivisions,
    "Number of times to subdivide the polygons in the surfaces." },
  { "-num_layers", ARGV_INT, (char*)0, (char*) &num_layers,
    "Number of intermediate layers to use." },
  { NULL, ARGV_END, NULL, NULL, NULL }
};

void measure_cortical_thickness( vector<polygons_struct> &layers,
				 ofstream &output) {
  Real total_dist, dist;
  int obj_index;
  Vector normal;
  Point origin_point, target_point;

  // David's method requires an object, so an object it'll get.
  object_struct *object = create_object( POLYGONS );

  // the surface we are starting from is the first one in the vector;
  vector<polygons_struct>::iterator origin = layers.begin();
  // and the next one is the first target		    
  vector<polygons_struct>::iterator target = layers.begin();
  target++;

  // now begin the search at each point of the first one
  unsigned int num_points = (*origin).n_points;
  for (unsigned int i=0; i < num_points; i++) {
    // set the point index to use. When casting from the first surface
    // it will correspond to i, later on it will correspond to the point
    // of intersection
    obj_index = i;
    total_dist = 0;
    // loop over all of the surfaces
    origin = layers.begin();
    target = layers.begin();
    target++;
    origin_point = (*origin).points[i];
    for ( ; target != layers.end(); target++ ) {
      object->specific.polygons = *target;

      // create a vector to cast the ray with
      /* t_{normal} segfaults for some reason
       *       SCALE_VECTOR( normal, (*origin).normals[obj_index], -1.0 );
       *       if( intersect_ray_with_object( &(*origin).points[obj_index],
       *                                      &normal, object,
       *                                      &obj_index, &dist, NULL ) == 0 ) {
       *         dist = 0.0;
       *       }
       */

      (void) find_closest_polygon_point( &origin_point, &(*target), 
					 &target_point );
      dist = distance_between_points( &origin_point, &target_point );
      origin_point = target_point;

      total_dist += dist;
      origin++;
    }
    output << total_dist << endl;
  }
    
}

int main(int argc, char *argv[]) {

  STRING surface1_filename, surface2_filename, thickness_filename;
  //int num_layers;
  int n_objects;
  File_formats format;
  polygons_struct *polygons1, *polygons2, *output;
  object_struct **object_list;
  object_struct *out_object;

  // argument handlingq
  if ( ParseArgv(&argc, argv, argTable, 0) || (argc != 4) ){
    cerr << "Usage: " << argv[0] << " [options] white_surface.obj grey_surface.obj thickness.txt" << endl;
    return(1);
  }
  surface1_filename  = argv[1];
  surface2_filename  = argv[2];
  thickness_filename = argv[3];

  // now open the surfaces
  if( input_graphics_file( surface1_filename, &format, &n_objects,
                           &object_list ) != OK || n_objects != 1 ||
      get_object_type(object_list[0]) != POLYGONS ) {
    cerr << "ERROR reading " << surface1_filename << endl;
  }
  polygons1 = get_polygons_ptr( object_list[0] );

  if( input_graphics_file( surface2_filename, &format, &n_objects,
                           &object_list ) != OK || n_objects != 1 ||
      get_object_type(object_list[0]) != POLYGONS ) {
    cerr << "ERROR reading " << surface2_filename << endl;
  }   
  polygons2 = get_polygons_ptr( object_list[0] );

  // open the output cortical thickness text file
  ofstream thickness_output(thickness_filename);
  if (! thickness_output) {
    cerr << "ERROR: could not open the output file " << thickness_filename
	 << endl;
    return(1);
  }

  // will hold the polygons for each of the layers
  vector<polygons_struct> layers(num_layers);

  // now create the layered surfaces
  if (num_layers > 0) {
    (void) create_intermediate_surfaces( polygons1, polygons2, 
					 num_layers, layers );
  }

  // add the second surface onto the vector so that it too can
  // be resampled
  layers.push_back(*polygons2);
  
  // subdivide the polygons if so desired
  if (num_subdivisions > 0) {
    (void) resample_polygon_vector(layers, num_subdivisions);
  }

  // create the bin trees
  create_polygon_bintrees( layers );

  // now add the starting surface - which will still have the original
  // number of polygons - to the vector
  layers.insert(layers.begin(), *polygons1);
  measure_cortical_thickness( layers, thickness_output );

//   out_object = create_object( POLYGONS );

//   // output the objects. Filenames will be of base_i.obj format
//   for (int i=0; i < num_layers; i++) {
//     char postfix[256];
//     STRING out_filename = create_string(output_filename_base);
//     sprintf(postfix, "_%d.obj", i);
//     out_filename = concat_strings( out_filename, postfix);
//     out_object->specific.polygons = layers[i];
//     output_graphics_file( out_filename, format, 1, &out_object );
//   }
    return 0;
}
