#include "laplacianGrid.h"

void main(int argc, char* argv[]) {
  vector<Real> xv, yv, zv;
  vector<Real>::iterator xit, yit, zit;


  laplacianGrid *grid = new laplacianGrid(argv[1], 0, 10000);

  //  grid->setVerbosity(10);

  cout << "Relaxing Equation." << endl;
  grid->relaxEquation(-1, 25);
  //  grid->output("grid.mnc");
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
  grid->computeAllThickness(0.1);

  grid->output(argv[2]);

}














