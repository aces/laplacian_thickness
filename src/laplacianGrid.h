#ifndef __LAPLACIANGRID__
#define __LAPLACIANGRID__

#include <mniLabelVolume.h>
#include <mniVolume.h>
#include <math.h>
#include <fstream>

using namespace std;

enum integrator { EULER, SECOND_ORDER_RK, FOURTH_ORDER_RK };
enum interpolation { NEAREST_NEIGHBOUR_INTERP = -1,
		     LINEAR_INTERP = 0,
		     CUBIC_INTERP = 2 };

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
  //! Holds the successive over-relaxation acceleration factor
  Real SOR;
  //! Holds a copy of the volume data in cache (for speed)
  Real * val;
  //! Hold the fixed grid
  mniVolume *fixedGrid;
  //! Hold the gradient volumes
  mniVolume *gradientX, *gradientY, *gradientZ;
  //! Hold the file whose values to average along streamline
  mniVolume *avgVolume;
  //! Variable to decide whether to compute average or length
  bool computeAverage;
  //! For progress reports
  progress_struct progressReport;
  //! Hold the level of verbosity
  int verbosity;
  //! Sets up the necessary volume information
  void initialiseVolumes(integrator type, nc_type gradientDataType);
  //! function pointer for the integration method to use.
  void (laplacianGrid::*integrationStep)(Real x, Real y, Real z,
					 Real dx, Real dy, Real dz, Real h,
                                         Real &newx,
                                         Real &newy,
                                         Real &newz );
  void eulerStep(Real x, Real y, Real z,
		 Real dx, Real dy, Real dz, Real h,
                 Real &newx,
                 Real &newy,
                 Real &newz );
  void fourthOrderRungeKuttaStep(Real x, Real y, Real z,
				 Real dx, Real dy, Real dz, Real h,
                                 Real &newx,
                                 Real &newy,
                                 Real &newz );
  void secondOrderRungeKuttaStep(Real x, Real y, Real z,
				 Real dx, Real dy, Real dz, Real h,
                                 Real &newx,
                                 Real &newy,
                                 Real &newz );

  void getDerivatives( Real x, Real y, Real z,
		       Real &dx, Real &dy, Real &dz );

  void setCacheVoxel( Real value, int v0, int v1, int v2 ) {
                      val[(v0*sizes[1]+v1)*sizes[2]+v2] = value; }
  Real getCacheVoxel( int v0, int v1, int v2 ) {
                      return( val[(v0*sizes[1]+v1)*sizes[2]+v2] ); }

  Real evaluate( Real x, Real y, Real z, interpolation interpType );  
  Real evaluateAvgVolume( Real x, Real y, Real z, interpolation interpType );  
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
		nc_type gradientDataType = NC_FLOAT);
  
  //! solve one iteration of laplace's equation
  float solveLaplace();
  //! relax the equation
  void relaxEquation(float convergenceCriteria, int maxIterations);
  //! create the normalised gradients volume in each direction
  void createNormalisedGradients();
  //! save the gradient volumes to file
  void saveGradients( char* gradientPrefix );
  //! compute average across streamline rather than length of streamline
  void setToAverageAlongStreamlines(char* avgFile,
				    nc_type volumeDataType=NC_SHORT);
  //! create a streamline at given voxel
  Real createStreamline(Real x0, Real y0, Real z0, Real h, 
			interpolation evalType);

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
