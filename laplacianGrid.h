#ifndef __LAPLACIANGRID__
#define __LAPLACIANGRID__

#include <mniLabelVolume.h>
#include <mniVolume.h>
#include "corticalMantle.h"

class laplacianGrid {
protected:
  //! Hold the fixed inner and outer boundaries
  int outerValue, innerValue;
  //! Hold the volume to solve over
  mniVolume volume;
  //! Hold the fixed grid
  mniLabelVolume fixedGrid;
  //! Hold the gradient volumes
  mniVolume gradientX, gradientY, gradientZ;
  
  //! code to initialise the gradient volumes
  void initialiseGradients();
public:
  //! constructor from corticalMantle
  laplacianGrid(corticalMantle *mantle);
  //! constructor from file
  laplacianGrid(char* mantleFile);
  
  //! solve one iteration of laplace's equation
  float solveLaplace();
  //! relax the mantle 
  void relaxMantle(float convergenceCriteria, int maxIterations);
  //! create a normalised gradient volume in each direction
  void createGradients();
  //! integrate gradients using using Euler's formula
  void createStreamlines();
  //! normalise according to units to get a thickness metric
  void createThicknessMetric();
  //! output the volume
  void output(char* filename);


};

#endif
