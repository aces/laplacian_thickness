#ifndef __LAPLACIANGRID__
#define __LAPLACIANGRID__

#include <mniLabelVolume.h>
#include <mniVolume.h>
#include <vector>
#include <math.h>
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
  //! For progress reports
  progress_struct progressReport;
  
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
  //! create a gradient volume in each direction
  void createGradients();
  //! normalise the gradient volumes
  void normaliseGradients();
  //! integrate gradients using Euler's formula
  void createStreamlines();
  //! normalise according to units to get a thickness metric
  void createThicknessMetric();
  //! create a streamline at given voxel
  void testFunction(vector<Real> &t);
  void createStreamline(int x0, int y0, int z0, int h, 
                        vector<Real> &Xvector,
                        vector<Real> &Yvector,
                        vector<Real> &Zvector);

  Real streamLength(vector<Real> &Xvector,
                    vector<Real> &Yvector,
                    vector<Real> &Zvector);

  void computeAllThickness();
  //  Real createStreamline(int x0, int y0, int z0, int h);
  //! output the volume
  void output(char* filename);



};

#endif
