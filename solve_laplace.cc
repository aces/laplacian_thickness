#include "solve_laplace.h"

void gradient_volume(Volume laplace_volume, int *sizes) {
  int v1, v2, v3;
  Real value_1, value_2, value_3;
  Volume gradient_1, gradient_2, gradient_3, tangent_field;

  // the gradient volume
  gradient_1 = copy_volume_definition(laplace_volume, NC_UNSPECIFIED,
					  FALSE, 0, 10000);
  /*
  gradient_2 = copy_volume_definition(laplace_volume, NC_UNSPECIFIED,
				      FALSE, 0, 10000);
  gradient_3 = copy_volume_definition(laplace_volume, NC_UNSPECIFIED,
				      FALSE, 0, 10000);
  */
  tangent_field = copy_volume_definition(laplace_volume, NC_UNSPECIFIED,
				      FALSE, 0, 10000);
  

  for (v1=1; v1 < sizes[0]-1; v1++) {
    for (v2=1; v2 < sizes[1]-1; v2++) {
      for (v3=1; v3 < sizes[2]-1; v3++) {
	if (get_volume_real_value(laplace_volume, v1, v2, v3, 0, 0) > 0) {
	  value_3 = 
	    (get_volume_real_value(laplace_volume, v1+1, v2, v3+1, 0, 0) -
	     get_volume_real_value(laplace_volume, v1-1, v2, v3-1, 0, 0)
	     / 2 );
	  value_2 = 
	    (get_volume_real_value(laplace_volume, v1, v2+1, v3, 0, 0) -
	     get_volume_real_value(laplace_volume, v1, v2-1, v3, 0, 0)
	     / 2 );
	  value_1 = 
	    (get_volume_real_value(laplace_volume, v1+1, v2, v3, 0, 0) -
	     get_volume_real_value(laplace_volume, v1-1, v2, v3, 0, 0)
	     / 2 );
	  set_volume_real_value(gradient_1, v1, v2, v3, 0, 0, value_1);

	  //	  set_volume_real_value(gradient_2, v1, v2, v3, 0, 0, value_2);
	  //	  set_volume_real_value(gradient_3, v1, v2, v3, 0, 0, value_3);
	  set_volume_real_value(tangent_field, v1, v2, v3, 0, 0,
				(sqrt(pow(value_1,2) + pow(value_2,2) +
				      pow(value_3,2))));
	}
      }
    }
  }

  delete_volume(laplace_volume);
  laplace_volume = copy_volume(tangent_field);
  delete_volume(gradient_1);
}
	  

/*
  Volume initialise_volume:

  Sets up a volume for subsequent solving. In other words it
  initialises the boundaries and gives a first rough guess for the
  intermediate pixels.  

  INPUT: 

  Volume label_volume: a label volume that has values for each of the
  surfaces as well as the intermediate mantle. This volume is modified
  in place with the new initialisation values.

*/

void initialise_volume(Volume label_volume) {
  int sizes[MAX_DIMENSIONS];
  int v1, v2, v3;
  Real value;

  get_volume_sizes(label_volume, sizes);
  
  for (v1=0; v1 < sizes[0]; v1++) {
    for (v2=0; v2 < sizes[1]; v2++) {
      for (v3=0; v3 < sizes[2]; v3++) {
	value = get_volume_real_value(label_volume, v1, v2, v3, 0, 0);
	switch((int)value) {
	case GREY_SURFACE:
	  //	  cout << "grey" << endl;
	  set_volume_real_value(label_volume, v1, v2, v3, 0, 0, OUTER_SURFACE);
	  break;
	case WHITE_SURFACE:
	  //	  cout << "white" << endl;
	  set_volume_real_value(label_volume, v1, v2, v3, 0, 0, INNER_SURFACE);
	  break;
	case BOTH_SURFACES:
	  //	  cout << "both" << endl;
	  set_volume_real_value(label_volume, v1, v2, v3, 0, 0, OUTER_SURFACE);
	  break;
	case MANTLE:
	  //	  cout << "mantle" << endl;
	  set_volume_real_value(label_volume, v1, v2, v3, 0, 0, MANTLE_INIT);
	  break;
	}
      }
    }
  }

}

void solve_laplace(Volume laplacian_volume, int *sizes) {
  int v1, v2, v3;
  Real initial_value, new_value;

  for (v1=1; v1 < sizes[0]-1; v1++) {
    for (v2=1; v2 < sizes[1]-1; v2++) {
      for (v3=1; v3 < sizes[2]-1; v3++) {
	initial_value = get_volume_real_value(laplacian_volume, v1, v2, 
					       v3, 0, 0);
	if (initial_value != INNER_SURFACE && 
	    initial_value != OUTER_SURFACE) {
	  new_value = (get_volume_real_value(laplacian_volume, v1+1, v2,
					     v3, 0, 0) + 
		       get_volume_real_value(laplacian_volume, v1, v2+1,
					     v3, 0, 0) + 
		       get_volume_real_value(laplacian_volume, v1, v2,
					     v3+1, 0, 0) + 
		       get_volume_real_value(laplacian_volume, v1-1, v2, 
					     v3, 0, 0) + 
		       get_volume_real_value(laplacian_volume, v1, v2-1,
					     v3, 0, 0) + 
		       get_volume_real_value(laplacian_volume, v1, v2,
					     v3-1, 0, 0) ) / 6;
	  set_volume_real_value(laplacian_volume, v1, v2, v3, 0, 0,
				new_value);
	}
      }
    }
  }
}
	  
