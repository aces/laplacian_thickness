#include "laplacianGrid.h"

extern "C" {
#include "ParseArgv.h"
}

using namespace std;

// set up argument parsing defaults
Real hValue         = 1;
int vValue          = 1;

// argument parsing table
ArgvInfo argTable[] = {
  { "-h", ARGV_FLOAT, (char *)0, (char *) &hValue,
    "H value to use for Eulerian integration" },
  { "-v", ARGV_INT, (char *)0, (char *) &vValue,
    "Verbosity level (higher = more verbose)" },

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
  grid->relaxEquation(-1, 25);
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
  grid->computeAllThickness( hValue );

  grid->output(argv[2]);
  return (0);

}














