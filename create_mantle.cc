#include <iostream.h>

extern "C" {
#include "volume_io.h"
#include "bicpl.h"
#include "ParseArgv.h"
}

#include "cortical_mantle.h"

// set up argument parsing stuff

int arg_relax_mantle = 1; // 1 for yes, 0 for no
int max_relax_iterations = 25;
int convergence_criteria = 4000;

ArgvInfo argTable[] = {
  { "-relax", ARGV_CONSTANT, (char*)1, (char*)&arg_relax_mantle,
    "find mantle through relaxation [default]" },
  { "-norelax", ARGV_CONSTANT, (char*)0, (char*)&arg_relax_mantle,
    "do not relax the mantle" },
  { "-max_relax_iterations", ARGV_INT, (char *)0, 
    (char *)&max_relax_iterations,
    "Maximum number of iterations for relaxation" },
  { "-convergence_criteria", ARGV_INT, (char *)0,
    (char *)&convergence_criteria,
    "Stop criteria for relaxation convergence" },
    
  { NULL, ARGV_END, NULL, NULL, NULL }
};


/********************
 void main
********************/
void main(int argc, char *argv[]) {
  Volume object_intersects, cortical_mantle, tmp_volume;
  char* dimnames[] = { MIzspace, MIxspace, MIyspace };
  int returnval, sizes[MAX_DIMENSIONS];
  jplPointMap point_map;


  if ( ParseArgv( &argc, argv, argTable, 0 ) || (argc != 3) ){
    cerr << "Usage: create_mantle [options] intersects.mnc out.mnc" << endl;
    exit(1);
  }

  char* intersects = argv[1];
  char* output = argv[2];


  // input the headers from the objects intersects file
  if (input_volume_header_only( intersects, 3, dimnames, &tmp_volume, 
				(minc_input_options *) NULL ) != OK ) {
    cerr << "ERROR: could not open object intersects" << endl;
    exit(1);
  }

  // open object intersects as a label volume
  object_intersects = create_label_volume( tmp_volume, NC_BYTE );

  if (load_label_volume(intersects, object_intersects) != OK) {
    cerr << "ERROR: could not open object intersects" << endl;
    exit(1);
  }

  get_volume_sizes(object_intersects, sizes);

  // find a few points within the mantle based on ray intersection
  find_cortical_mantle(object_intersects);

  // find the remaining points in the mantle through neigbour checks
  // and relaxation. This step can be turned off from the command line.
  if (arg_relax_mantle == 1) {
    cout << "now relaxing mantle" << endl;
    relax_mantle(object_intersects, sizes, max_relax_iterations,
		 convergence_criteria);
  }

  cout << "Now writing map to volume" << endl;
  save_label_volume(output, intersects, object_intersects, 10);

  delete_volume(object_intersects);



}

