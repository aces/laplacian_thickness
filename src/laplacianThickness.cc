#include "laplacianGrid.h"

extern "C" {
#include "ParseArgv.h"
#include "createLaplacianGrid.h"
}

using namespace std;

enum operation { FROM_GRID, FROM_SURFACES };

// set up argument parsing defaults
Real hValue                      = 0.1;        /* in voxel space */
Real convergence                 = 0.00001;
int  maxIterations               = 300;
int  vValue                      = 1;          // verbosity level
char *objFile                    = NULL;
int  potentialOnly               = 0;
integrator integration           = EULER;      // EULER, SECOND_ORDER_RK, FOURTH_ORDER_RK
operation mode                   = FROM_SURFACES;
int  include_white_boundary      = 0;
int  include_grey_boundary       = 0;
char *likeFile                   = NULL; 
char *averageFile                = NULL;
nc_type volumeType               = NC_SHORT;   // NC_BYTE, NC_SHORT, NC_INT, NC_FLOAT, NC_DOUBLE
nc_type gradientsType            = NC_FLOAT;   // MUST have float or double (CL)
interpolation boundaryType       = NEAREST_NEIGHBOUR_INTERP;   // MUST be nearest neighbour (CL)
                                   // NEAREST_NEIGHBOUR_INTERP, LINEAR_INTERP, CUBIC_INTERP

nc_type gradientsJunk            = NC_FLOAT;                   // for backward compatibility (CL)
interpolation boundaryJunk       = NEAREST_NEIGHBOUR_INTERP;   // for backward compatibility (CL)

// argument parsing table
ArgvInfo argTable[] = {
  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nOperation Type:" },
  { "-from_surfaces", ARGV_CONSTANT, (char *)FROM_SURFACES,
    (char *) &mode,
    "Create the grid from two surfaces [Default]" },
  { "-from_grid", ARGV_CONSTANT, (char *)FROM_GRID, (char *) &mode,
    "Use a prepared grid." },

  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nCortical Mantle Options:" },
  { "-include_white_boundary", ARGV_CONSTANT, (char *)1, 
    (char *) &include_white_boundary,
    "Force the inclusion of the white-matter surface  \n\t\t\t\tboundary. [Default: false]" },
  { "-include_grey_boundary", ARGV_CONSTANT, (char *)1,
    (char *) &include_grey_boundary,
    "Force the inclusion of the grey-matter surface  \n\t\t\t\tboundary. [Default: false]"},
  { "-like", ARGV_STRING, (char *)0, (char *) &likeFile,
    "MINC file whose shape the cortical mantle will take" },

  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nLaplacian Equation solving options:" },
  { "-h", ARGV_FLOAT, (char *)0, (char *) &hValue,
    "H value to use for Eulerian integration" },
  { "-v", ARGV_INT, (char *)0, (char *) &vValue,
    "Verbosity level (higher = more verbose)" },
  { "-convergence", ARGV_FLOAT, (char *)0, (char *) &convergence,
    "Stop criteria for equation relaxation." },
  { "-max_iterations", ARGV_INT, (char *)0, (char *) &maxIterations,
    "Maximum number of iterations for relaxation." },
  { "-object_eval", ARGV_STRING, (char *)0, (char *) &objFile,
    "Evaluate thickness only at vertices of the obj file. \n\t\t\t\tOutput text rather than minc." },
  { "-potential_only", ARGV_CONSTANT, (char *)1, (char *) &potentialOnly,
    "Output only the potential field and stop." },
  { "-average_along_streamlines", ARGV_STRING, (char *)0, 
    (char *) &averageFile,
    "Compute mean value of voxels in specified file along streamlines." },

  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nIntegration Options:" },
  { "-euler", ARGV_CONSTANT, (char *)EULER, (char *)&integration,
    "Use Euler method for integration (Default)" },
  { "-2nd_rk", ARGV_CONSTANT, (char *)SECOND_ORDER_RK, (char *)&integration,
    "Use second order Runge-Kutta integration" },
  { "-4th_rk", ARGV_CONSTANT, (char *)FOURTH_ORDER_RK, (char *)&integration,
    "Use fourth order Runge-Kutta integration" },

  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nData-type control:" },
  { "-volume-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&volumeType,
    "Volume is byte format." },
  { "-volume-short", ARGV_CONSTANT, (char *)NC_SHORT, (char *)&volumeType,
    "Volume is short integer format." },
  { "-volume-int", ARGV_CONSTANT, (char *)NC_INT, (char *)&volumeType,
    "Volume is integer format." },
  { "-volume-float", ARGV_CONSTANT, (char *)NC_FLOAT, (char *)&volumeType,
    "Volume is floating point precision." },
  { "-volume-double", ARGV_CONSTANT, (char *)NC_DOUBLE, (char *)&volumeType,
    "Volume is double floating point precision." },

  { "-gradients-byte", ARGV_CONSTANT, (char *)NC_BYTE, (char *)&gradientsJunk,
    "Gradients in byte format are no longer supported. Using float by default." },
  { "-gradients-short", ARGV_CONSTANT, (char *)NC_SHORT,
    (char *)&gradientsJunk,
    "Gradients in short integer format are no longer supported. Using float by default." },
  { "-gradients-int", ARGV_CONSTANT, (char *)NC_INT, (char *)&gradientsJunk,
    "Gradients in integer format are no longer supported. Using float by default." },
  { "-gradients-float", ARGV_CONSTANT, (char *)NC_FLOAT,
    (char *)&gradientsType,
    "Gradients are floating point precision." },
  { "-gradients-double", ARGV_CONSTANT, (char *)NC_DOUBLE,
    (char *)&gradientsType,
    "Gradients are double floating point precision." },

  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "\nBoundary evaluation control:" },
  { "-boundary-nearest-neighbour", ARGV_CONSTANT,
    (char *)NEAREST_NEIGHBOUR_INTERP,
    (char *)&boundaryType,
    "Use nearest neighbour interpolation." },
  { "-boundary-linear", ARGV_CONSTANT, (char *)LINEAR_INTERP,
    (char *)&boundaryJunk,
    "Linear interpolation is no longer supported. Using nearest neighbour by default." },
  { "-boundary-cubic", ARGV_CONSTANT, (char *)CUBIC_INTERP,
    (char *)&boundaryJunk,
    "Cubic interpolation is no longer supported. Using nearest neighbour by default." },

  { NULL, ARGV_END, NULL, NULL, NULL }
};

int main(int argc, char* argv[]) {

  if ( ParseArgv(&argc, argv, argTable, 0) ) {
    cerr << endl << "Usage:\t" << argv[0] << " [options] -like sample.mnc grey_surface.obj \n\twhite_surface.obj output_thickness.mnc" << endl;
    cerr << "\tor" << endl;
    cerr << "Usage:\t" << argv[0] << " [options] -like sample.mnc -object_eval surface.obj\n\t grey_surface.obj white_surface.obj output_thickness.txt" << endl << endl;
    return (1);
  }

  // After parsing, check new value or argc
  if( ( mode == FROM_SURFACES && argc != 4 ) ||
      ( mode == FROM_GRID && argc != 3 ) ) {
    cerr << endl << "Usage:\t" << argv[0] << " [options] -like sample.mnc grey_surface.obj \n\twhite_surface.obj output_thickness.mnc" << endl;
    cerr << "\tor" << endl;
    cerr << "Usage:\t" << argv[0] << " [options] -like sample.mnc -object_eval surface.obj\n\t grey_surface.obj white_surface.obj output_thickness.txt" << endl << endl;
    return (1);
  }

  int outside_value = 10;
  int inside_value = 0;
  char *grey_surface;
  char *white_surface;
  char *out_filename;
  char *grid_file;
  laplacianGrid *grid;
   

  if (mode == FROM_SURFACES) {

    grey_surface = argv[1];
    white_surface = argv[2];
    out_filename = argv[3];
    
    cout << "Creating cortical mantle." << endl;
    Volume input_grid;
    if( create_mantle( likeFile, grey_surface, white_surface,
                       include_white_boundary, include_grey_boundary,
                       &input_grid ) != OK ) {
      cerr << "Error: Cannot create mantle" << endl << endl;
      return (1);
    } else {
      cout << "Mantle created" << endl;
    }
    grid = new laplacianGrid(input_grid, 
                             inside_value,
                             outside_value,
                             integration,
			     volumeType,
			     gradientsType);

  }
  else if (mode == FROM_GRID) {
    grid_file = argv[1];
    out_filename = argv[2];
    grid = new laplacianGrid(grid_file, 
                             inside_value,
                             outside_value,
                             integration,
			     volumeType,
			     gradientsType);
    cout << "after constructor " << out_filename << endl;

  }


  grid->setVerbosity( vValue );

  cout << "Relaxing Equation." << endl;
  grid->relaxEquation(convergence, maxIterations);

  if ( potentialOnly == 1 ) {
    grid->output(out_filename);
    exit(0);
  }

  cout << "Creating normalised gradients." << endl;
  grid->createNormalisedGradients();

  if ( averageFile == NULL ) {
    cout << "Beginning computation of thicknesses." << endl;
  }
  else {
    grid->setToAverageAlongStreamlines( averageFile );
    cout << "Beginning averaging of values along streamlines." << endl;
  }

  

  if( objFile == NULL ) {
    grid->computeAllThickness( hValue, boundaryType );
    grid->output(out_filename);
  }
  else {
    cout << "  Only computing at each vertex" << endl;
    grid->computeAllThickness( hValue, objFile, boundaryType );
    grid->output( out_filename, true );
  }

  return (0);

}

