extern "C" {
#include "volume_io.h"
#include "bicpl.h"
}

#include <iostream.h>
#include "solve_laplace.h"

void main(int argc, char *argv[]) {
  Volume volume;
  minc_output_options options;
  int sizes[MAX_DIMENSIONS];
  int i = 0;
  

  if ( input_volume(argv[1], 3, NULL, NC_FLOAT, FALSE, 0.0, 10000.0, TRUE,
		    &volume, (minc_input_options *)NULL) != OK ) {
    cerr << "ERROR: could not open volume " << argv[1] << endl;
    exit(1);
  }

  initialise_volume(volume);

  get_volume_sizes(volume, sizes);
  
  while (i < 15) {
    solve_laplace(volume, sizes);
    i++;
  }
  gradient_volume(volume, sizes);

  set_default_minc_output_options(&options);
  set_minc_output_real_range(&options, 0.0, 10000.0);

  output_volume("laplace_gradient.mnc", NC_FLOAT, FALSE, 0.0, 10000.0, volume,
		"laplacian test volume", &options);
			 
}

