#include "laplacianGrid.h"

extern "C" {
#include "ParseArgv.h"
#include "createLaplacianGrid.h"
}

using namespace std;

// set up argument parsing defaults
Real hValue                      = 0.1;
Real convergence                 = 0.00001;
int  maxIterations               = 300;
int  vValue                      = 1;
char *objFile                    = NULL;
int  potentialOnly               = 0;
int  include_white_boundary      = 0;
int  include_grey_boundary       = 0;
char *likeFile                   = NULL; 

// argument parsing table
ArgvInfo argTable[] = {
  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "Cortical Mantle Options:" },
  { "-include_white_boundary", ARGV_CONSTANT, (char *)1, 
    (char *) &include_white_boundary,
    "Force the inclusion of the white-matter surface boundary. [Default: false]" },
  { "-include_grey_boundary", ARGV_CONSTANT, (char *)1,
    (char *) &include_grey_boundary,
    "Force the inclusion of the grey-matter surface boundary. [Default: false]"},
  { "-like", ARGV_STRING, (char *)0, (char *) &likeFile,
    "MINC file whose shape the cortical mantle will take" },

  { NULL, ARGV_HELP, (char *)NULL, (char *)NULL,
    "Laplacian Equation solving options:" },
  { "-h", ARGV_FLOAT, (char *)0, (char *) &hValue,
    "H value to use for Eulerian integration" },
  { "-v", ARGV_INT, (char *)0, (char *) &vValue,
    "Verbosity level (higher = more verbose)" },
  { "-convergence", ARGV_FLOAT, (char *)0, (char *) &convergence,
    "Stop criteria for equation relaxation." },
  { "-max_iterations", ARGV_INT, (char *)0, (char *) &maxIterations,
    "Maximum number of iterations for relaxation." },
  { "-object_eval", ARGV_STRING, (char *)0, (char *) &objFile,
    "Evaluate thickness only at vertices of the obj file. Output text rather than minc." },
  { "-potential_only", ARGV_CONSTANT, (char *)1, (char *) &potentialOnly,
    "Output only the potential field and stop (Default: False)" },

  { NULL, ARGV_END, NULL, NULL, NULL }
};

int main(int argc, char* argv[]) {
  vector<Real> xv, yv, zv;
  vector<Real>::iterator xit, yit, zit;

  if ( ParseArgv(&argc, argv, argTable, 0) || (argc != 4 )) {
    cerr << "Usage: " << argv[0] << " [options] -like sample.mnc grey_surface.obj white_surface.obj output_thickness.mnc" << endl;
    cerr << "or" << endl;
    cerr << "Usage: " << argv[0] << " [options] -like sample.mnc -object_eval surface.obj grey_surface.obj white_surface.obj output_thickness.txt" << endl;

    return (1);
  }

  char *grey_surface = argv[1];
  char *white_surface = argv[2];
  char *out_filename = argv[3];

  int outside_value = 10000;
  int inside_value = 0;

  cout << "Creating cortical mantle." << endl;
  Volume input_grid = create_mantle( likeFile,
				     grey_surface,
				     white_surface,
				     include_white_boundary,
				     include_grey_boundary );
  cout << "Mantle created" << endl;

  laplacianGrid *grid = new laplacianGrid(input_grid, 
					  inside_value,
					  outside_value);

  grid->setVerbosity( vValue );

  cout << "Relaxing Equation." << endl;
  grid->relaxEquation(convergence, maxIterations);

  if ( potentialOnly == 1 ) {
    grid->output(out_filename);
    exit(0);
  }

  cout << "Creating gradients." << endl;
  grid->createGradients();
  cout << "Normalising gradients." << endl;
  grid->normaliseGradients();

  /*
  cout << "test" << endl;
    grid->createStreamline(123,110,35,1,xv,yv,zv);
    
    xit = xv.begin();
    yit = yv.begin();
    zit = zv.begin();
    while (xit != xv.end()) {
    cout << *xit << " " << *yit << " " << *zit << endl;
    xit++; yit++; zit++;
    }
    cout << "Size: " << xv.size() << endl;
    cout << "Length: " << grid->streamLength(xv, yv, zv) << endl;
  */
  //  grid->setVerbosity(10);
  //  grid->createStreamline(94,85,154, 0.5, xv, yv, zv);
  //  grid->createStreamline(83,172,91, 1, xv, yv, zv);

  cout << "Beginning computation of thicknesses." << endl;

  if( objFile == NULL ) {
    grid->computeAllThickness( hValue );
    grid->output(out_filename);
  }
  else {
    cout << "  Only computing at each vertex" << endl;
    grid->computeAllThickness( hValue, objFile );
    grid->output( out_filename, true );
  }

  return (0);

}














