#include "cortical_mantle.h"

#include "cortical_mantle.h"


/*
  check_if_between_surfaces - tests whether a voxel lies between the inner 
  and outer surface. The strategy is to expand a line out in each of the
  cardinal directions to see which surface is first intersected. If it
  intersects both white and grey surfaces it is considered to be in
  the cortical mantle.

  NOTE: this method will not find all points within the mantle. It's goal
  instead is to find a point known to lie within the mantle and to then
  expand out from that point to fill in the rest of the cortex.  */

int check_if_between_surfaces(Volume object_intersects,
			      int* sizes, 
			      int coord1, int coord2, int coord3, 
			      int longest_line) {

  int grey_intersect, white_intersect, long_line;
  int value;
  grey_intersect = 0;
  white_intersect = 0;

  // build an array out of the coordinates to make looping easier
  int coords[3] = { coord1, coord2, coord3 };

  // cast a line into all the cardinal directions, stopping when you
  // hit a surface or the end of the volume

  for (int i=0; i<=2; i++) {
    // first in descending order, i.e. coords[i]--
    long_line = coords[i] - longest_line;
    while (coords[i] >= long_line && coords[i] != 0) {
      value = get_volume_label_data_5d(object_intersects, coords[0], 
					 coords[1], coords[2], 0, 0);
      if (value == BOTH_SURFACES) { 
	break;
      }
      else if (value == GREY_SURFACE) {
	grey_intersect++;
	break;
      }
      else if (value == WHITE_SURFACE) {
	white_intersect++;
	break;
      }
      coords[i]--;
    }
    //reinitialise the coordinates to their default values
    coords[0] = coord1;
    coords[1] = coord2;
    coords[2] = coord3;

    long_line = coords[i] + longest_line;
    //and now in ascending order, i.e. coords[i]++
    while (coords[i] <= long_line && coords[i] != sizes[i]) {
      //    cout << "coord: " << coords[i] <<endl;
      value = get_volume_label_data_5d(object_intersects, coords[0], 
					 coords[1], coords[2], 0, 0);
      //      if (value > 1)
      //	cout << "value in intersects: " << value << endl;
      
      if (value == BOTH_SURFACES) {
	break;
      }
      else if (value == GREY_SURFACE) {
	grey_intersect++;
	break;
      }
      else if (value == WHITE_SURFACE) {
	white_intersect++;
	//	cout << "white: " << white_intersect << endl;
	break;
      }
      coords[i]++;
    }
    //reinitialise the coordinates to their default values
    coords[0] = coord1;
    coords[1] = coord2;
    coords[2] = coord3;
  }

  // determine the return value based on the number of grey and white
  // intersects
  //  if (white_intersect > 0 && grey_intersect > 0)
  //   cout << "White: " << white_intersect << 
  //  " Grey: "  << grey_intersect << endl;

  if (grey_intersect > 2 && white_intersect > 2) {
    //cout << "returned 1" << endl;
    return 1;
  }
  else {
    //    cout << "returned 0" << endl;
    return 0;
  }
}

/*
  find_cortical_mantle: 
  takes two arguments: a label volume, which will be modified, and a regular
  volume which contains the surface volume intersect data. The label volume
  will have four possible values at the end: 
  0 for no interest
  1 for grey surface
  2 for white surface
  3 for point in the cortical mantle
*/

int find_cortical_mantle(Volume object_intersects) {
  int sizes[MAX_DIMENSIONS];
  int v1, v2, v3, is_in_mantle;
  int value;
  progress_struct progress;


  // run through the volume and determine the mantle
  get_volume_sizes(object_intersects, sizes);
  initialize_progress_report(&progress, FALSE, sizes[0], 
			     "Computing cortical mantle");
  for (v1 = 0; v1 < sizes[0]; v1++) {
    update_progress_report(&progress, v1 + 1);
    for (v2 = 0; v2 < sizes[1]; v2++) {
      for (v3 = 0; v3 < sizes[2]; v3++) {
	value = get_volume_label_data_5d(object_intersects, v1, v2, v3, 0, 0);

	if (value == NO_VALUE) { 
	  if ( check_if_between_surfaces(object_intersects, sizes, 
					 v1, v2, v3, 5) == 1) {
	    set_volume_label_data_5d(object_intersects, v1, v2, v3, 
				     0, 0, MANTLE);
	  }
        }
      }
    }
  }
  terminate_progress_report(&progress);

  return 1; //ok, so there's no error checking yet ...
}
