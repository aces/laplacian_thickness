#include "surface_functions.h"

void resample_polygon_vector( vector<polygons_struct> &layers,
                              unsigned short subsampling_steps ) {

  // in order to maintain tetrahedral topology, the new number of items
  // must be the old times a multiple of four.
  int i;
  vector<polygons_struct>::iterator j;
  for (i=0; i< subsampling_steps; i++) {
    unsigned short new_n_items = layers[0].n_items * 4;
    for (j = layers.begin();j != layers.end(); j++ ) {
      subdivide_polygons( &*j );
    }
  }
  for (j = layers.begin(); j != layers.end(); j++) {
    compute_polygon_normals( &*j );
  }
}

/*! Creates the polygon bintrees for all the surfaces
 */
void create_polygon_bintrees( vector<polygons_struct> &layers ) {
  vector<polygons_struct>::iterator j;
  for (j = layers.begin(); j != layers.end(); j++) {
    create_polygons_bintree( &*j, ROUND( (Real) (*j).n_items * 0.4 ));
  }
}

/*! Creates a weighted averages from X surfaces
 */
void create_intermediate_surfaces(polygons_struct *surface1,
                                  polygons_struct *surface2,
                                  unsigned short num_layers,
                                  vector<polygons_struct> &layers) {
  
  // will hold the weights for each averaging run
  unsigned short weight1, weight2;

  // ensure that averaging is legal
  if (! polygons_are_same_topology( surface1, surface2 ) ) {
    cerr << "ERROR: polygons have to be of same topology" << endl
         << "Error in create_intermediate_surface" << endl;
  }
  
  // loop over each of the layers
  for (unsigned short n=0; n < num_layers; n++) {
    // create the output polygon
    (void) copy_polygons( surface1, &layers[n] );

    // determine the weights for this layer
    weight1 = 1 + n;
    weight2 = num_layers - n;
      
    // do the weighted averaging
    for (int i=0; i< surface1->n_points; i++) {
      // extract the individual points
      Point s1 = surface1->points[i];
      Point s2 = surface2->points[i];
      Point is;

      // compute the weighted avg of each point
      Point_x(s1) *= weight1;
      Point_y(s1) *= weight1;
      Point_z(s1) *= weight1;
      
      Point_x(s2) *= weight2;
      Point_y(s2) *= weight2;
      Point_z(s2) *= weight2;

      Point_x(is) = ( Point_x(s1) + Point_x(s2) ) / ( weight1 + weight2 );
      Point_y(is) = ( Point_y(s1) + Point_y(s2) ) / ( weight1 + weight2 );
      Point_z(is) = ( Point_z(s1) + Point_z(s2) ) / ( weight1 + weight2 );
      
      //      cout << "Before point assignment" << endl;
      layers[n].points[i] = is;
    }
    compute_polygon_normals( &layers[n] );
  }
}
