#ifndef __LAPLACIAN2DGRID__
#define __LAPLACIAN2DGRID__

#include <mniLabelVolume.h>
#include <mniVolume.h>
#include <vector>
#include <algo.h>
#include <math.h>
#include <fstream>
#include "laplacianGlobals.h"

using namespace std;

/*
enum integrator { EULER, SECOND_ORDER_RK, FOURTH_ORDER_RK };
enum interpolation { NEAREST_NEIGHBOUR_INTERP = -1,
		     LINEAR_INTERP = 0,
		     CUBIC_INTERP = 2 };
*/
static STRING XYdimOrder[] = {MIxspace, MIyspace};

//! a class to hold the two dimensional index
class twodindex {
 public:
  int x;
  int y;
  twodindex(int x,int y) {this->x = x; this->y = y;};
};

class laplacian2DGrid {
protected:
  //! Hold the volume sizes
  int *sizes;
  //! Hold the fixed inner and outer boundaries
  int outerValue, innerValue;
  //! Hold the volume to solve over
  mniVolume *volume;
  //! The text file to be output to if solving only at vertices
  FILE *outputVertexFile;
  //! Holds the thickness info at each vertex
  Real *thicknessPerVertex;
  //! Holds number of vertices in obj file
  int numVertices;
  //! Hold the fixed grid
  mniLabelVolume *fixedGrid;
  //! Hold the gradient volumes
  mniVolume *gradientX, *gradientY;
  //! For progress reports
  progress_struct progressReport;
  //! Hold the level of verbosity
  int verbosity;
  //! Sets up the necessary volume information
  void initialiseVolumes(integrator type, nc_type gradientDataType);
  //! function pointer for the integration method to use.
  void (laplacian2DGrid::*integrationStep)(vector<Real> &Xvector,
                                         vector<Real> &Yvector,
                                         Real dx, Real dy, Real h,
                                         Real &newx,
                                         Real &newy,
                                         int currentStep);
  void eulerStep(vector<Real> &Xvector,
                 vector<Real> &Yvector,
                 Real dx, Real dy, Real h,
                 Real &newx,
                 Real &newy,
                 int currentStep);
  void getDerivatives( Real x, Real y, Real &dx, Real &dy );
  Real evaluate( Real x, Real y, interpolation interpType );  

public:

  //! constructor from file
  laplacian2DGrid(char* mantleFile,
                int innerValue,
                int outerValue,
		integrator type = EULER,
		nc_type volumeDataType = NC_SHORT,
		nc_type gradientDataType = NC_BYTE);
  //! constructor from volume_struct pointer
  laplacian2DGrid(Volume mantleVolume,
		int innerValue,
		int outerValue,
		integrator type = EULER,
		nc_type volumeDataType = NC_SHORT,
		nc_type gradientDataType = NC_BYTE);
 //! solve one iteration of laplace's equation
  float solveLaplace();
  //! relax the equation
  void relaxEquation(float convergenceCriteria, int maxIterations);
  //! create a gradient volume in each direction
  void createGradients();
  //! an initial guess for the solving based on distances.
  void initialGridGuess(int sparsity=5);
  //! normalise the gradient volumes
  void normaliseGradients();
  //! integrate gradients using Euler's formula
  void createStreamlines();
  //! normalise according to units to get a thickness metric
  void createThicknessMetric();
  //! create a streamline at given voxel
  void testFunction(vector<Real> &t);
  void createStreamline(Real x0, Real y0, Real h, 
                        vector<Real> &Xvector,
                        vector<Real> &Yvector,
                        interpolation evalType);
  Real streamLength(vector<Real> &Xvector,
                    vector<Real> &Yvector);
  void computeAllThickness(Real h, 
			   interpolation evalType = NEAREST_NEIGHBOUR_INTERP);
  //! output the volume
  void output(char* filename, bool isTextFile=false);
  //! set the level of verbosity
  void setVerbosity(int newVerbosity) { this->verbosity = newVerbosity; };


};

#endif
