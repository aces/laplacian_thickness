#ifndef __LAPLACIANGRID__
#define __LAPLACIANGRID__

#include <mniLabelVolume.h>
#include <mniVolume.h>
#include "corticalMantle.h"

class laplacianGrid {
protected:
  //! Hold the volume sizes
  int *sizes;
  //! Hold the fixed inner and outer boundaries
  int outerValue, innerValue;
  //! Hold the volume to solve over
  mniVolume *volume;
  //! Hold the fixed grid
  mniLabelVolume *fixedGrid;
  //! Hold the gradient volumes
  mniVolume *gradientX, *gradientY, *gradientZ;
  
public:
  //! constructor from corticalMantle
  laplacianGrid(corticalMantle *mantle);
  //! constructor from file
  laplacianGrid(char* mantleFile,
                int innerValue,
                int outerValue);
  
  //! solve one iteration of laplace's equation
  float solveLaplace();
  //! relax the equation
  void relaxEquation(float convergenceCriteria, int maxIterations);
  //! create a normalised gradient volume in each direction
  void createGradients();
  //! integrate gradients using using Euler's formula
  void createStreamlines();
  //! normalise according to units to get a thickness metric
  void createThicknessMetric();
  //! create a streamline at given voxel
  Real createStreamline(int x0, int y0, int z0, int h);
  //! output the volume
  void output(char* filename);


};

#endif
