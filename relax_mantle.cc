#include "relax_mantle.h"


/**********************************
  determine_neighbouring_mantle_points:
  
  The assumption is that a few points in the mantle have already been
  found through other means (such as ray intersection). The task of
  this routine is therefore to find each of the points in the mantle
  and to check whether its neighbours also deserve to be in the
  mantle. The check is very simple: if the neighbour has no value
  assigned to it, it is in the mantle - otherwise it would have the
  value of one of the two surfaces.
  
  Inputs:
     Volume cortical_mantle: the volume containing distinct values for
     the two surfaces and those points already having been identified
     as lying in the cortex

     int *sizes: an array of three values containing the sizes of the
     cortical_mantle volume

  Returns the number of neighbouring points that were changed - this
  value can be used to test for relaxation.
**********************************/

int determine_neighbouring_mantle_points(Volume cortical_mantle, int *sizes) {
  int indices[3], tmp_indices[3];
  int i, j, num_values_changed;

  num_values_changed = 0;

  // loop over the volume, checking for neighbours
  for (indices[0]=1; indices[0] < sizes[0]-1; indices[0]++) {
    for (indices[1]=1; indices[1] < sizes[1]-1; indices[1]++) {
      for (indices[2]=1; indices[2] < sizes[2]-1; indices[2]++) {
	// check if point is in mantle
	if ( get_volume_label_data_5d(cortical_mantle, indices[0], indices[1],
				      indices[2], 0, 0) == MANTLE ) {
	  // check all neighbours - if value is 0 change value to 4
	  for (i=0; i < 3; i++) {
	    for (j=-1; j < 2; j++) {
	      tmp_indices = indices;
	      tmp_indices[i] += j;
	      if ( get_volume_label_data_5d(cortical_mantle, tmp_indices[0],
					    tmp_indices[1], tmp_indices[2],
					    0, 0) == NO_VALUE ) {
		/* assing a temporary number to neighbouring points
                   within the mantle. The only reason a temporary
                   value is used is that otherwise the image will
                   quickly be considered to be entirely mantle if any
                   of the original points were incorrectly labeled as
                   mantle points. This way the damage will be more
                   localised */
 		set_volume_label_data_5d(cortical_mantle, tmp_indices[0],
					 tmp_indices[1], tmp_indices[2],
					 0, 0, TMP_MANTLE);
 		num_values_changed++;
	      }
	    }
	  }
	}
      }
    }
  }

  /* the code here replaces the temporary value with the correct value
     for a mantle point. See above for an explanation of this lunacy */
  for (indices[0]=0; indices[0] < sizes[0]; indices[0]++) {
    for (indices[1]=0; indices[1] < sizes[1]; indices[1]++) {
      for (indices[2]=0; indices[2] < sizes[2]; indices[2]++) {
	if ( get_volume_label_data_5d(cortical_mantle, indices[0],
				      indices[1], indices[2],
				      0, 0) == TMP_MANTLE ) {
	  set_volume_label_data_5d(cortical_mantle, indices[0],
				   indices[1], indices[2],
				   0, 0, MANTLE);
	}
      }
    }
  }

  cout << "Num changed: " << num_values_changed << endl;
  return num_values_changed;
}

/******************************************
  relax_mantle:
  
  This routine calls determine_neighbouring_mantle_points repeatedly
  until the system relaxes itself or the maximum number of allowed
  iterations is reached.

  Inputs:

    Volume cortical_mantle: the volume containing the surface and
    mantle points

    int *sizes: an array of three values containing the surface sizes

    int max_iterations: the maximum number of iterations allowed
    before giving up on relaxation

    int convergence_criteria: the algorithm will stop if the number of
    voxels changed dips below this number
***************************************/

void relax_mantle(Volume cortical_mantle, int *sizes, int max_iterations,
		  int convergence_criteria) {
  int i=0;
  int last_iterations[3] = {0, 0, 0};
  int result;
  while (i < max_iterations) {
    result = determine_neighbouring_mantle_points(cortical_mantle, sizes);
    if (result < convergence_criteria)
      break;
    
    // push iterations onto stack of last iterations
    // not yet implemented

    i++;
  }
}
