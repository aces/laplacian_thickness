#include "laplacianGrid.h"

void main(int argc, char* argv[]) {
  vector<Real> xv, yv, zv;
  vector<Real>::iterator xit, yit, zit;


  laplacianGrid *grid = new laplacianGrid(argv[1], 0, 10000);

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


  cout << "Beginning computation of thicknesses." << endl;
  grid->computeAllThickness();

  grid->output(argv[2]);

}














