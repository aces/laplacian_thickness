#include "laplacianGrid.h"

void main() {
  laplacianGrid *grid = new laplacianGrid("/home/jason/data/sendai/idac_aoba_00934_mantle.mnc", 0, 10000);
  grid->relaxEquation(-1, 15);
  grid->createGradients();
  grid->output("testgrid.mnc");
}
