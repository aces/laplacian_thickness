#include "laplacianGrid.h"

void main(int argc, char* argv[]) {
  vector<Real> xv(3), yv(3), zv(3);

  laplacianGrid *grid = new laplacianGrid(argv[1], 0, 10000);
  grid->relaxEquation(-1, 15);
  grid->createGradients();
  grid->normaliseGradients();
  grid->createStreamline(96,108,123,1,xv,yv,zv);
  cout << xv[1] << endl;
  grid->output("testgrid.mnc");
}
