/***
 * A revised t_normal thickness implementation
 * $Id$
 ***/

extern "C" {
#include <bicpl.h>
#include <ParseArgv.h>
}

#include "surface_functions.h"

#include <iostream>
#include <fstream>

using namespace std;

// argument parsing stuff

Real layered_normals_thickness( vector <polygons_struct> &layers,
				int index ) {
  Real distance;
  int new_index;
  vector<Vector> all_normals;
  vector<polygons_struct>::iterator j, source, target;
  Vector normals;

  int num_layers = 0;
  for(j = layers.begin(); j != layers.end(); j++) {
    SCALE_VECTOR( normals, (*j).normals[index], -1.0 );
    all_normals.push_back(normals);
    num_layers++;
  }
  vector<Vector>::iterator i;
  Vector_x(normals) = 0;
  Vector_y(normals) = 0;
  Vector_z(normals) = 0;
  for (i = all_normals.begin(); i != all_normals.end(); i++) {
    Vector_x(normals) += Vector_x(*i);
    Vector_y(normals) += Vector_y(*i);
    Vector_z(normals) += Vector_z(*i);
  }
  Vector_x(normals) = Vector_x(normals) / num_layers;
  Vector_y(normals) = Vector_y(normals) / num_layers;
  Vector_z(normals) = Vector_z(normals) / num_layers;

  object_struct *object = create_object( POLYGONS );
  source = layers.begin();
  target = layers.end();
  target--;
  object->specific.polygons = *target;
  if (intersect_ray_with_object( &(*source).points[index],
				 &normals, object, &new_index,
				 &distance, NULL ) == 0 ) {
    distance = 0.0;
  }
  return distance;
}

Real measure_normal_thickness( vector<polygons_struct> &layers,
			       int index ) {
  Real distance;
  int new_index;
  vector<Real> all_distances;
  Vector normals;
  vector<polygons_struct>::iterator j, finish, target;
  finish = layers.end();
  finish--;
  for (j = layers.begin(); j != finish; j++) {
    target = j;
    target++;
    object_struct *object = create_object( POLYGONS );
    object->specific.polygons = *target;

    SCALE_VECTOR( normals, (*j).normals[index], -1.0 );
    if (intersect_ray_with_object( &(*j).points[index],
				   &normals, object, &new_index,
				   &distance, NULL ) == 0 ) {
      distance = 0.0;
    }
    all_distances.push_back(distance);
    //    cout << "dist: " << distance << endl;
  }
  vector<Real>::iterator i;
  for (i = all_distances.begin(); i != all_distances.end(); i++ ) {
    distance += pow(*i, 2);
  }
  distance = sqrt(distance);
  return distance;
}

int main(int argc, char *argv[]) {
  STRING surface1_filename, surface2_filename, thickness_filename;
  int n_objects;
  File_formats format;
  polygons_struct *polygons1, *polygons2;
  object_struct **object_list;
  
  surface1_filename = argv[1];
  surface2_filename = argv[2];
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
  int num_layers = 8;

  vector<polygons_struct> layers(num_layers);

  create_intermediate_surfaces( polygons1, polygons2, 
				num_layers, layers );
  layers.push_back( *polygons2 );
  layers.insert( layers.begin(), *polygons1 );
  create_polygon_bintrees( layers );
  for (int i=0; i < polygons1->n_points; i++) {
    //    Real dist = measure_normal_thickness( layers, i );
    Real dist = layered_normals_thickness( layers, i );
    
    thickness_output << dist << endl;
  }
  thickness_output.close();
}
