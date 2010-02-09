/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
#include "laplacianGrid.h"
#include <string>

// constructor from file
laplacianGrid::laplacianGrid(char* mantleFile, 
                             int innerValue, 
                             int outerValue,
                             integrator integrationType,
			     nc_type volumeDataType,
			     nc_type gradientDataType) {


  this->innerValue = innerValue;
  this->outerValue = outerValue;

  if( getenv( "LAPLACE_SOR" ) ) {
    sscanf( getenv("LAPLACE_SOR"), "%lf", &(this->SOR) );
    if( this->SOR <= 0.0 || this->SOR >= 2.0 ) this->SOR = 1.75;
  } else {
    this->SOR = 1.75;   // optimal value for 3-D Laplacian (Claude)
  }
  this->val = NULL;

  int ndims;
  STRING * dimension_names;
  get_file_dimension_names( mantleFile, &ndims, &dimension_names );
  if( ndims != 3 ) {
    cerr << "Error: Dimension of input volume " << mantleFile <<
            " must be 3." << endl;
    exit( 1 );
  }

  // create the grid volume from file
  this->fixedGrid = new mniVolume(mantleFile, 0.0, 0.0, 3, dimension_names);

  // construct the volume from the mantle file, but signed
  this->volume = new mniVolume(mantleFile,
			       0.0,
			       0.0,
                               3,
                               dimension_names,
                               volumeDataType,
                               TRUE,
                               TRUE,
                               NULL);

  // initialise the sizes member
  this->sizes = new int[3];
  this->sizes = this->fixedGrid->getSizes();

  this->volume->setRealRange( (Real)this->innerValue, (Real)this->outerValue );

  this->initialiseVolumes(integrationType, gradientDataType);
}

// constructor from volume_struct
laplacianGrid::laplacianGrid(Volume mantleVolume,
			     int innerValue,
			     int outerValue,
			     integrator integrationType,
			     nc_type volumeDataType,
			     nc_type gradientDataType) {
  this->innerValue = innerValue;
  this->outerValue = outerValue;

  if( getenv( "LAPLACE_SOR" ) ) {
    sscanf( getenv("LAPLACE_SOR"), "%lf", &(this->SOR) );
    if( this->SOR <= 0.0 || this->SOR >= 2.0 ) this->SOR = 1.75;
  } else {
    this->SOR = 1.75;   // optimal value for 3-D Laplacian (Claude)
  }
  this->val = NULL;

  // take the volume_struct pointer as the grid
  this->fixedGrid = new mniVolume(mantleVolume);
  //this->fixedGrid->output("grid.mnc");
  this->volume = new mniVolume(this->fixedGrid, TRUE, volumeDataType, TRUE);
  this->volume->setRealRange( (Real)this->innerValue, (Real)this->outerValue );

  // initialise the sizes member
  this->sizes = new int[3];
  this->sizes = this->fixedGrid->getSizes();

  // volume definition copy won't copy values, while a direct copy
  // will keep the original data-type. So this is the ugly compromise.
  for (int v0=0; v0 < this->sizes[0]; v0++) {
    for (int v1=0; v1 < this->sizes[1]; v1++) {
      for (int v2=0; v2 < this->sizes[2]; v2++) {
        this->volume->setVoxel( this->fixedGrid->getVoxel(v0, v1, v2),
                                v0, v1, v2);
      }
    }
  }

  this->initialiseVolumes(integrationType, gradientDataType);
}
    

void laplacianGrid::initialiseVolumes(integrator integrationType,
				      nc_type gradientDataType) {

  // and construct the gradient volumes - using volume definition copy
  this->gradientX = new mniVolume(this->volume, TRUE, gradientDataType, TRUE);
  this->gradientY = new mniVolume(this->volume, TRUE, gradientDataType, TRUE);
  this->gradientZ = new mniVolume(this->volume, TRUE, gradientDataType, TRUE);


  // set their real ranges - safe here since I don't care about
  // the input volume
  Real lower = -1.0;
  Real upper = 1.0;
  gradientX->setRealRange(lower, upper);
  gradientY->setRealRange(lower, upper);
  gradientZ->setRealRange(lower, upper);

  // assume that length of streamlines will be computed
  this->computeAverage = false;

  // set default verbosity to 0
  this->verbosity = 0;
  
  // set the function pointer
  cout << "Using integration type: ";
  if (integrationType == EULER) {
    this->integrationStep = &laplacianGrid::eulerStep;
    cout << "Euler." << endl;
  }
  else if (integrationType == SECOND_ORDER_RK) {
    this->integrationStep = &laplacianGrid::secondOrderRungeKuttaStep;
    cout << "second order Runge-Kutta." << endl;
  }
  else if (integrationType == FOURTH_ORDER_RK) {
    this->integrationStep = &laplacianGrid::fourthOrderRungeKuttaStep;
    cout << "fourth order Runge-Kutta." << endl;
  }

}

// (in voxel space)
inline void laplacianGrid::getDerivatives( Real x, Real y, Real z,
					   Real &dx, Real &dy, Real &dz ) {

  dx = this->gradientX->getInterpolatedVoxel(x,y,z,CUBIC_INTERP);
  dy = this->gradientY->getInterpolatedVoxel(x,y,z,CUBIC_INTERP);
  dz = this->gradientZ->getInterpolatedVoxel(x,y,z,CUBIC_INTERP);
}

// (in voxel space)
Real laplacianGrid::evaluate( Real x, Real y, Real z, 
			     interpolation interpType) {
  //  return this->fixedGrid->getInterpolatedVoxel(x,y,z);
  return this->fixedGrid->getInterpolatedVoxel(x, y, z, interpType);
}

Real laplacianGrid::evaluateAvgVolume( Real x, Real y, Real z,
				       interpolation interpType ) {
  return this->avgVolume->getInterpolatedVoxel(x, y, z, interpType);
}
								   
void laplacianGrid::setToAverageAlongStreamlines(char* avgFile, 
						 nc_type volumeDataType) {
  this->computeAverage = true;

  int ndims;
  STRING * dimension_names;
  get_file_dimension_names( avgFile, &ndims, &dimension_names );
  if( ndims != 3 ) {
    cerr << "Error: Dimension of input volume " << avgFile <<
            " must be 3." << endl;
    exit( 1 );
  }

  this->avgVolume = new mniVolume(avgFile,
				  0.0,
				  0.0,
				  3,
                                  dimension_names,
				  volumeDataType,
				  TRUE,
				  TRUE,
				  NULL);
}



// take an integration step using Euler's method (in voxel space)
inline void laplacianGrid::eulerStep(Real x, Real y, Real z,
                                     Real dx, Real dy, Real dz, Real h, 
                                     Real &newx,
                                     Real &newy,
                                     Real &newz ) {

  newx = x + dx * h;
  newy = y + dy * h;
  newz = z + dz * h;

}

// take an integration step using a second order Runge-Kutta Model (in voxel space)
inline void laplacianGrid::secondOrderRungeKuttaStep(
                                          Real x, Real y, Real z,
                                          Real dx, Real dy, Real dz, Real h, 
                                          Real &newx,
                                          Real &newy,
                                          Real &newz ) {

  Real hh = h * 0.5;
  Real kx, ky, kz;
  Real ddx, ddy, ddz;

  kx = x + dx * hh;
  ky = y + dy * hh;
  kz = z + dz * hh;

  this->getDerivatives(kx, ky, kz, ddx, ddy, ddz);

  newx = x + ddx * h;
  newy = y + ddy * h;
  newz = z + ddz * h;
}
  

// take an integration step using a fourth order Runge-Kutta Model (in voxel space)
inline void laplacianGrid::fourthOrderRungeKuttaStep(
					  Real x, Real y, Real z,
					  Real dx, Real dy, Real dz, Real h, 
                                          Real &newx,
                                          Real &newy,
                                          Real &newz ) {

  Real ddx[3], ddy[3], ddz[3];
  Real kx, ky, kz;

  Real hh = h*0.5;

  // take first step
  kx = x + dx * hh;
  ky = y + dy * hh;
  kz = z + dz * hh;

  // second step
  this->getDerivatives(kx, ky, kz, ddx[0], ddy[0], ddz[0]);

  kx = x + ddx[0] * hh;
  ky = y + ddy[0] * hh;
  kz = z + ddz[0] * hh;
  
  // third step
  this->getDerivatives(kx, ky, kz, ddx[1], ddy[1], ddz[1]);
  
  // fourth step
  kx = x + ddx[1] * h;
  ky = y + ddy[1] * h;
  kz = z + ddz[1] * h;

  this->getDerivatives(kx, ky, kz, ddx[2], ddy[2], ddz[2]);
  
  newx = x + h * (dx + 2.0*(ddx[0]+ddx[1]) + ddx[2]) / 6.0;
  newy = y + h * (dy + 2.0*(ddy[0]+ddy[1]) + ddy[2]) / 6.0;
  newz = z + h * (dz + 2.0*(ddz[0]+ddz[1]) + ddz[2]) / 6.0;

}


// one iteration of solving laplace's equation
float laplacianGrid::solveLaplace() {

  Real residual = 0.0;

  Real initialValue, newValue, oldValue, tmpValue;

  for (int v0=1; v0 < this->sizes[0]-1; v0++) {
    for (int v1=1; v1 < this->sizes[1]-1; v1++) {
      for (int v2=1; v2 < this->sizes[2]-1; v2++) {
        initialValue = this->fixedGrid->getVoxel(v0,v1,v2);
        if (initialValue > (this->innerValue +0.1) &&
            initialValue < (this->outerValue -0.1) ) {
          tmpValue = (this->getCacheVoxel(v0+1, v1, v2) + 
                      this->getCacheVoxel(v0, v1+1, v2) +
                      this->getCacheVoxel(v0, v1, v2+1) +
                      this->getCacheVoxel(v0-1, v1, v2) +
                      this->getCacheVoxel(v0, v1-1, v2) +
                      this->getCacheVoxel(v0, v1, v2-1))/6.0;

          // Apply SOR to speed up convergence
          oldValue = this->getCacheVoxel(v0, v1, v2);
          newValue = oldValue + SOR * ( tmpValue - oldValue );

          // Disable SOR if newValue out of bounds (SOR=1).
          if( newValue < this->innerValue ) {
            newValue = tmpValue;
          } else if( newValue > this->outerValue ) {
            newValue = tmpValue;
          }

          residual += ABS( oldValue - newValue );  // residual

          this->setCacheVoxel(newValue, v0, v1, v2);

        }
      }
    }
  }

  return residual;
}
        
// relax the equation
void laplacianGrid::relaxEquation(float convergenceCriteria,
                                  int maxIterations) {

  float convergence, normalizeFactor;

  int currentIteration = 0;

  // Copy the data in cache for speed
  this->val = new Real[this->sizes[0]*this->sizes[1]*this->sizes[2]];
  for (int v0=0; v0 < this->sizes[0]; v0++) {
    for (int v1=0; v1 < this->sizes[1]; v1++) {
      for (int v2=0; v2 < this->sizes[2]; v2++) {
        this->setCacheVoxel( this->volume->getVoxel(v0, v1, v2), v0, v1, v2 );
      }
    }
  }

  do {
    currentIteration++;
    convergence = this->solveLaplace();

    if( currentIteration == 1 ) normalizeFactor = 1.0 / convergence;
    convergence *= normalizeFactor;  // normalize by initial residual

    if (this->verbosity >= 1) {
      cout << "Converging: iteration # " << currentIteration 
           << " with convergence result of " << convergence << endl;
    }
  } while( convergence > convergenceCriteria &&
           currentIteration < maxIterations );
  
  // Copy back the data from cache to volume after solving
  for (int v0=0; v0 < this->sizes[0]; v0++) {
    for (int v1=0; v1 < this->sizes[1]; v1++) {
      for (int v2=0; v2 < this->sizes[2]; v2++) {
        this->volume->setVoxel( this->getCacheVoxel( v0, v1, v2 ), v0, v1, v2 );
      }
    }
  }

  delete [] this->val;
  this->val = NULL;
}

// create a normalised gradient volume in each of the cardinal directions
//       normalise the X gradient at that position like so:
//       Nx = dx / [dx^2 + dy^2 + dz^2]^0.5
//       by dividing with the vector magnitude.
//
void laplacianGrid::createNormalisedGradients() {

  Real nx, ny, nz, dx, dy, dz, mag;

  // create the gradients using central difference formula,
  // e.g. dx = volume(x+1,y,z) - volume(x-1,y,z,) / 2
  for (int x=1; x < this->sizes[0]-1; x++) {
    for (int y=1; y < this->sizes[1]-1; y++) {
      for (int z=1; z < this->sizes[2]-1; z++) {
        if (this->fixedGrid->getVoxel(x,y,z) > this->innerValue &&
            this->fixedGrid->getVoxel(x,y,z) < this->outerValue ) {

          // note: no need to divide by 2 since we normalise
          dx = this->volume->getVoxel(x+1,y,z) - this->volume->getVoxel(x-1,y,z);
          dy = this->volume->getVoxel(x,y+1,z) - this->volume->getVoxel(x,y-1,z);
          dz = this->volume->getVoxel(x,y,z+1) - this->volume->getVoxel(x,y,z-1);
          mag = sqrt( dx*dx + dy*dy + dz*dz );
          if( mag > 1.0e-06 ) {
            nx = dx / mag;
            ny = dy / mag;
            nz = dz / mag;
          } else {
            nx = 0.0;
            ny = 0.0;
            nz = 0.0;
          }

          this->gradientX->setVoxel( nx, x,y,z );
          this->gradientY->setVoxel( ny, x,y,z );
          this->gradientZ->setVoxel( nz, x,y,z );
          
        } else {
          this->gradientX->setVoxel(0, x,y,z);
          this->gradientY->setVoxel(0, x,y,z);
          this->gradientZ->setVoxel(0, x,y,z);
        }
      }
    }
  }
}

void laplacianGrid::saveGradients(char* gradientPrefix) {
  string gradientXFilename;
  string gradientYFilename;
  string gradientZFilename;

  gradientXFilename = gradientPrefix;
  gradientXFilename.append("_dx.mnc");
  gradientYFilename = gradientPrefix;
  gradientYFilename.append("_dy.mnc");
  gradientZFilename = gradientPrefix;
  gradientZFilename.append("_dz.mnc");

  cout << "Outputting X gradient to " << gradientXFilename << endl;
  this->gradientX->output((char *)gradientXFilename.c_str());

  cout << "Outputting Y gradient to " << gradientYFilename << endl;
  this->gradientY->output((char *)gradientYFilename.c_str());

  cout << "Outputting Z gradient to " << gradientZFilename << endl;
  this->gradientZ->output((char *)gradientZFilename.c_str());

}

// use integration method, first towards one then towards the other surface
// Return the length of the streamline (in voxel units).
// Alternately returns the average values along a streamline. This behaviour
// is enabled by calling setToAverageAlongStreamlines.
Real laplacianGrid::createStreamline(Real x0, Real y0, Real z0, Real h,
				     interpolation evalType) {

  Real oldx, oldy, oldz, newx, newy, newz;
  Real dx, dy, dz, mag, length, sumValue;
  int nsteps;

  length = 0;
  sumValue = 0;
  nsteps = 0;

  Real eval0 = this->evaluate( x0, y0, z0, evalType );
  Real evaluation = eval0;

  // move towards outside surface first
  oldx = x0;
  oldy = y0;
  oldz = z0;
  while (evaluation < this->outerValue) {

    this->getDerivatives( oldx, oldy, oldz, dx, dy, dz );
    mag = dx*dx + dy*dy + dz*dz;
    // test for stagnation point (point not moving due to zero derivative)
    if( mag < 1.0e-6 ) {
      evaluation = this->outerValue;
    } else {
      (*this.*integrationStep)( oldx, oldy, oldz, dx, dy, dz, h,
                                newx, newy, newz );

      if (computeAverage) {
	sumValue += this->evaluateAvgVolume( newx, newy, newz, evalType );
	nsteps++;
      }
      length += sqrt( (newx-oldx)*(newx-oldx) + (newy-oldy)*(newy-oldy) +
                      (newz-oldz)*(newz-oldz) );
      
      // find out where we are
      evaluation = this->evaluate( newx, newy, newz, evalType );

      // This is a quick check to make sure that the streamline
      // is not bifurcating in saddle points or getting trapped
      // around a local max (isolated csf voxel inside gray
      // matter). If the streamline length is more than 4 times
      // the straight line distance between the starting point
      // and the end point, then likely something is wrong.

      Real line_dist = sqrt( (newx-x0)*(newx-x0) + (newy-y0)*(newy-y0) +
                             (newz-z0)*(newz-z0) );

      if( length > 4.0 * line_dist ) {
        evaluation = this->outerValue;
      }

      oldx = newx;
      oldy = newy;
      oldz = newz;
    }
  }

  // move towards inner surface
  oldx = x0;
  oldy = y0;
  oldz = z0;
  evaluation = eval0;
  Real length1 = length;
  length = 0.0;

  while (evaluation > this->innerValue) {

    this->getDerivatives( oldx, oldy, oldz, dx, dy, dz );
    mag = dx*dx + dy*dy + dz*dz;
    // test for stagnation point (point not moving due to zero derivative)
    if( mag < 1.0e-6 ) {
      evaluation = this->innerValue;
    } else {
      (*this.*integrationStep)( oldx, oldy, oldz, dx, dy, dz, -h,
                                newx, newy, newz );

      if (computeAverage) {
	sumValue += this->evaluateAvgVolume( newx, newy, newz, evalType );
	nsteps++;
      }
      length += sqrt( (newx-oldx)*(newx-oldx) + (newy-oldy)*(newy-oldy) +
		      (newz-oldz)*(newz-oldz) );

      // find out where we are
      evaluation = this->evaluate( newx, newy, newz, evalType );

      // This is a quick check to make sure that the streamline
      // is not bifurcating in saddle points or getting trapped
      // around a local max (isolated csf voxel inside gray
      // matter). If the streamline length is more than 4 times
      // the straight line distance between the starting point
      // and the end point, then likely something is wrong.

      Real line_dist = sqrt( (newx-x0)*(newx-x0) + (newy-y0)*(newy-y0) +
                             (newz-z0)*(newz-z0) );

      if( length > 4.0 * line_dist ) {
        evaluation = this->innerValue;
      }

      oldx = newx;
      oldy = newy;
      oldz = newz;
    }
  }
  length += length1;

  if (computeAverage) {
    //cout << sumValue << " " << nsteps<< " " << sumValue / nsteps << endl;

    return (sumValue / nsteps );
  }
  else {
    // return( length );
    return( length1 / length );
  }
}

void laplacianGrid::computeAllThickness( Real h, char *objFile,
					 interpolation evalType ) {
  File_formats 	format;
  int 		point, nPoints, nObjects;
  Point 	*points;
  object_struct **objects;
  Real          *voxel, result;
  Real          interpVoxel[3];
		
  // open the file
  if( input_graphics_file( objFile, &format, &nObjects, &objects )
      != OK ) {
    cerr << "ERROR: could not open file " << objFile << endl;
    exit(1);
  }
  
  if( nObjects != 1 ) {
    cerr << "WARNING: more than one object in " << objFile << endl;
  }
  
  nPoints = get_object_points( objects[0], &points );
  this->thicknessPerVertex = new Real[nPoints];
  this->numVertices = nPoints;

  // normalise unit to be millimetres.
  // the step size ought to be uniform. If not we're screwed anyway.
  Real separations[MAX_DIMENSIONS];
  get_volume_separations( this->fixedGrid->getVolume(), separations );
  
  initialize_progress_report( &this->progressReport, FALSE, nPoints,
			      "Computing streamlines" );
  for (int i=0; i < nPoints; i++) {
    voxel = this->fixedGrid->convertWorldToVoxel(
						 RPoint_x(points[i]),
						 RPoint_y(points[i]),
						 RPoint_z(points[i]) );
    
    if (this->fixedGrid->getInterpolatedVoxel(voxel[0], voxel[1], voxel[2]) 
	< this->outerValue &&
	this->fixedGrid->getInterpolatedVoxel
	(voxel[0], 
	 voxel[1], 
	 voxel[2]) 
	> this->innerValue) {
      result = this->createStreamline(voxel[0], voxel[1], voxel[2], h, evalType );
      if (computeAverage == false){
	result *= separations[0];
      }
    } else {
      result = 0;
      if ( this->verbosity > 1 ) {
        cout << "Could not evaluate at vertex " << i << endl;
      }
    }
    this->thicknessPerVertex[i] = result;
    update_progress_report( &this->progressReport, i );
  }
}
	

void laplacianGrid::computeAllThickness(Real h, interpolation evalType) {

  // normalise unit to be millimetres.
  // the step size ought to be uniform. If not we're screwed anyway.
  Real separations[MAX_DIMENSIONS];
  get_volume_separations( this->fixedGrid->getVolume(), separations );
  
  initialize_progress_report(&this->progressReport, FALSE, this->sizes[0]-1,
                             "Computing thickness streamlines");

  for (int v1=1; v1 < this->sizes[0]-1; v1++) {
    for (int v2=1; v2 < this->sizes[1]-1; v2++) {
      for (int v3=1; v3 < this->sizes[2]-1; v3++) {
        if (this->fixedGrid->getVoxel(v1, v2, v3) > this->innerValue && 
            this->fixedGrid->getVoxel(v1, v2, v3) < this->outerValue) {
          Real result = this->createStreamline(v1, v2, v3, h, evalType);
          result = ( 1.0 - result ) * this->outerValue;
	  if (computeAverage == false) {
	    result *= separations[0];
	  }
          this->volume->setVoxel(result ,v1, v2, v3);
        }
        else {
          this->volume->setVoxel(0, v1, v2, v3);
        }
      }
    }
    update_progress_report(&this->progressReport, v1);
  }
  terminate_progress_report(&this->progressReport);
}

void laplacianGrid::output(char *filename, bool isTextFile) {

  if ( isTextFile == false ) {
    this->volume->output(filename);
  }
  else {

    // create the output file to write the vertex info to
    ofstream outfile(filename);
    for ( int i=0; i < this->numVertices; ++i ) {
      outfile << this->thicknessPerVertex[i] << endl;
    }
    outfile.close();
  }
}
