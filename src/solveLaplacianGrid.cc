#include "laplacianGrid.h"

extern "C" {
#include "ParseArgv.h"
}

using namespace std;

// set up argument parsing defaults
Real hValue         = 0.1;
Real convergence    = 0.00001;
int maxIterations   = 300;
int vValue          = 1;
char *objFile       = NULL;

// argument parsing table
ArgvInfo argTable[] = {
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

  { NULL, ARGV_END, NULL, NULL, NULL }
};

int main(int argc, char* argv[]) {
  vector<Real> xv, yv, zv;
  vector<Real>::iterator xit, yit, zit;

  if ( ParseArgv(&argc, argv, argTable, 0) || (argc != 3 )) {
    cerr << "Usage: " << argv[0] << " [options] input_grid.mnc output_thickness.mnc" << endl;
    return (1);
  }

  laplacianGrid *grid = new laplacianGrid(argv[1], 0, 10000);

  grid->setVerbosity( vValue );

  cout << "Relaxing Equation." << endl;
  grid->relaxEquation(convergence, maxIterations);
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
    grid->output(argv[2]);
  }
  else {
    cout << "  Only computing at each vertex" << endl;
    grid->computeAllThickness( hValue, objFile );
    grid->output( argv[2], true );
  }

  return (0);

}














