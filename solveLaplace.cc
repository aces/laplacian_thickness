#include "laplacianGrid.h"

void main(int argc, char* argv[]) {
  vector<Real> xv, yv, zv;
  vector<Real>::iterator xit, yit, zit;


  laplacianGrid *grid = new laplacianGrid(argv[1], 0, 10000);
  grid->relaxEquation(-1, 25);
  grid->output("grid.mnc");
  grid->createGradients();
  grid->normaliseGradients();
  grid->createStreamline(129,119,22,1,xv,yv,zv);

  xit = xv.begin();
  yit = yv.begin();
  zit = zv.begin();
  while (xit != xv.end()) {
    cout << *xit << " " << *yit << " " << *zit << endl;
    xit++; yit++; zit++;
  }
  cout << "Size: " << xv.size() << endl;
  cout << "Length: " << grid->streamLength(xv, yv, zv) << endl;

  grid->computeAllThickness();

  grid->output("thickness.mnc");

}
