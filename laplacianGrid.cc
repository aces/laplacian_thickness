#include "laplacianGrid.h"

// initialise the two volumes used for solving

// constructor from file
laplacianGrid::laplacianGrid(char* mantleFile) {

  // NOTE: uses the default Z-X-Y dimension ordering
  this->fixedGrid = new mniLabelVolume(mantleFile);
