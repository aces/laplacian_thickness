#ifndef __LAPLACIANGRID__
#define __LAPLACIANGRID__

#include <mniLabelVolume.h>
#include <mniVolume.h>
#include "laplacianGlobals.h"
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

class laplacianGrid {
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
  mniVolume *gradientX, *gradientY, *gradientZ;
  //! For progress reports
  progress_struct progressReport;
  //! Hold the level of verbosity
  int verbosity;
  //! Sets up the necessary volume information
  void initialiseVolumes(integrator type, nc_type gradientDataType);
  //! function pointer for the integration method to use.
  void (laplacianGrid::*integrationStep)(vector<Real> &Xvector,
					 vector<Real> &Yvector,
					 vector<Real> &Zvector,
					 Real dx, Real dy, Real dz, Real h,
                                         Real &newx,
                                         Real &newy,
                                         Real &newz,
					 int currentStep);
  void eulerStep(vector<Real> &Xvector,
		 vector<Real> &Yvector,
		 vector<Real> &Zvector,
		 Real dx, Real dy, Real dz, Real h,
                 Real &newx,
                 Real &newy,
                 Real &newz,
		 int currentStep);
  void fourthOrderRungeKuttaStep(vector<Real> &Xvector,
				 vector<Real> &Yvector,
				 vector<Real> &Zvector,
				 Real dx, Real dy, Real dz, Real h,
                                 Real &newx,
                                 Real &newy,
                                 Real &newz,
				 int currentStep);
  void secondOrderRungeKuttaStep(vector<Real> &Xvector,
				 vector<Real> &Yvector,
				 vector<Real> &Zvector,
				 Real dx, Real dy, Real dz, Real h,
                                 Real &newx,
                                 Real &newy,
                                 Real &newz,
				 int currentStep);

  void getDerivatives( Real x, Real y, Real z,
		       Real &dx, Real &dy, Real &dz );

  Real evaluate( Real x, Real y, Real z, interpolation interpType );  
public:
  //! constructor from corticalMantle
  //  laplacianGrid(corticalMantle *mantle);
  //! constructor from file
  laplacianGrid(char* mantleFile,
                int innerValue,
                int outerValue,
		integrator type = EULER,
		nc_type volumeDataType = NC_SHORT,
		nc_type gradientDataType = NC_BYTE);
  //! constructor from volume_struct pointer
  laplacianGrid(Volume mantleVolume,
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
  //! normalise the gradient volumes
  void normaliseGradients();
  //! integrate gradients using Euler's formula
  void createStreamlines();
  //! normalise according to units to get a thickness metric
  void createThicknessMetric();
  //! create a streamline at given voxel
  void testFunction(vector<Real> &t);
  void createStreamline(Real x0, Real y0, Real z0, Real h, 
                        vector<Real> &Xvector,
                        vector<Real> &Yvector,
                        vector<Real> &Zvector,
			interpolation evalType);

  Real streamLength(vector<Real> &Xvector,
                    vector<Real> &Yvector,
                    vector<Real> &Zvector);

  void computeAllThickness(Real h, 
			   interpolation evalType = NEAREST_NEIGHBOUR_INTERP);
  //! overloaded to only create streamlines at vertices
  void computeAllThickness(Real h, char* objFile,
			   interpolation evalType = NEAREST_NEIGHBOUR_INTERP);
  //  Real createStreamline(int x0, int y0, int z0, int h);
  //! output the volume
  void output(char* filename, bool isTextFile=false);
  //! set the level of verbosity
  void setVerbosity(int newVerbosity) { this->verbosity = newVerbosity; };


};

#endif
