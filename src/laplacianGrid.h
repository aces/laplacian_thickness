#ifndef __LAPLACIANGRID__
#define __LAPLACIANGRID__

#include <mniLabelVolume.h>
#include <mniVolume.h>
#include <vector>
#include <math.h>
#include <fstream>
//#include "corticalMantle.h"

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
  mniVolume *fixedGrid;
  //! Hold the gradient volumes
  mniVolume *gradientX, *gradientY, *gradientZ;
  //! For progress reports
  progress_struct progressReport;
  //! Hold the level of verbosity
  int verbosity;
  //! Sets up the necessary volume information
  void initialiseVolumes();
  //! function pointer for the integration method to use.
  void (laplacianGrid::*integrationStep)(vector<Real> &Xvector,
			  vector<Real> &Yvector,
			  vector<Real> &Zvector,
			  Real dx, Real dy, Real dz, Real h,
			  vector<Real>::iterator XinsertIt,
			  vector<Real>::iterator YinsertIt,
			  vector<Real>::iterator ZinsertIt,
			  int currentStep);
  void eulerStep(vector<Real> &Xvector,
			  vector<Real> &Yvector,
			  vector<Real> &Zvector,
			  Real dx, Real dy, Real dz, Real h,
			  vector<Real>::iterator XinsertIt,
			  vector<Real>::iterator YinsertIt,
			  vector<Real>::iterator ZinsertIt,
			  int currentStep);
  void rungeKuttaStep(vector<Real> &Xvector,
			  vector<Real> &Yvector,
			  vector<Real> &Zvector,
			  Real dx, Real dy, Real dz, Real h,
			  vector<Real>::iterator XinsertIt,
			  vector<Real>::iterator YinsertIt,
			  vector<Real>::iterator ZinsertIt,
			  int currentStep);
  void secondOrderRungeKuttaStep(vector<Real> &Xvector,
			  vector<Real> &Yvector,
			  vector<Real> &Zvector,
			  Real dx, Real dy, Real dz, Real h,
			  vector<Real>::iterator XinsertIt,
			  vector<Real>::iterator YinsertIt,
			  vector<Real>::iterator ZinsertIt,
			  int currentStep);


  void getDerivatives( Real x, Real y, Real z,
		       Real &dx, Real &dy, Real &dz );

  
public:
  //! constructor from corticalMantle
  //  laplacianGrid(corticalMantle *mantle);
  //! constructor from file
  laplacianGrid(char* mantleFile,
                int innerValue,
                int outerValue);
  //! constructor from volume_struct pointer
  laplacianGrid(Volume mantleVolume,
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
  void createStreamline(int x0, int y0, int z0, Real h, 
                        vector<Real> &Xvector,
                        vector<Real> &Yvector,
                        vector<Real> &Zvector);

  Real streamLength(vector<Real> &Xvector,
                    vector<Real> &Yvector,
                    vector<Real> &Zvector);

  void computeAllThickness(Real h);
  //! overloaded to only create streamlines at vertices
  void computeAllThickness(Real h, char* objFile);
  //  Real createStreamline(int x0, int y0, int z0, int h);
  //! output the volume
  void output(char* filename, bool isTextFile=false);
  //! set the level of verbosity
  void setVerbosity(int newVerbosity) { this->verbosity = newVerbosity; };


};

#endif
