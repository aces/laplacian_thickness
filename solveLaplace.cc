#include "laplacianGrid.h"

void main(int argc, char* argv[]) {
  laplacianGrid *grid = new laplacianGrid(argv[1], 0, 10000);
  grid->relaxEquation(-1, 15);
  grid->createGradients();
  grid->normaliseGradients();
  grid->output("testgrid.mnc");
}
