#include "laplacian2DGrid.h"

// constructor from file
laplacian2DGrid::laplacian2DGrid( char *gridFile,
                                  int innerValue,
                                  int outerValue,
                                  integrator integrationType,
                                  nc_type volumeDataType,
                                  nc_type gradientDataType ) {
  
  this->innerValue = innerValue;
  this->outerValue = outerValue;

  // create the grid volume from file
  fixedGrid = new mniLabelVolume(gridFile, 0.0, 0.0, 3, XYZdimOrder);
  // construct the volume from the mantle file, but signed
  volume = new mniVolume( gridFile, 0.0, 0.0, 3, XYZdimOrder,
                          volumeDataType, TRUE, TRUE, NULL );

  // initialise the sizes
  sizes = new int[MAX_DIMENSIONS];
  sizes = fixedGrid->getSizes();

  volume->setRealRange(0, 20);
  initialiseVolumes(integrationType, gradientDataType);

}

void laplacian2DGrid::initialiseVolumes(integrator integrationType,
                                        nc_type gradientDataType) {
  
  // construct the gradient volumes using volume definition copy
  gradientX = new mniVolume(volume, TRUE, gradientDataType, TRUE);
  gradientY = new mniVolume(volume, TRUE, gradientDataType, TRUE);
  
  // set real range
  Real lower = -10;
  Real upper = 10;
  gradientX->setRealRange(lower, upper);
  gradientY->setRealRange(lower, upper);

  // set default verbosity to 0;
  verbosity = 0;

  // set the function pointer - note, no actual selection like in the
  // 3D case happens yet - that still has to be written.
  integrationStep = &laplacian2DGrid::eulerStep;
}

inline void laplacian2DGrid::getDerivatives( Real x, Real y,
                                             Real &dx, Real &dy ) {
  dx = gradientX->getInterpolatedVoxel(x,y,0,2);
  dy = gradientY->getInterpolatedVoxel(x,y,0,2);
}

Real laplacian2DGrid::evaluate(Real x, Real y, interpolation interpType) {
  return fixedGrid->getInterpolatedVoxel(x,y,0,interpType);
}

inline void laplacian2DGrid::eulerStep(vector<Real> &Xvector,
                                       vector<Real> &Yvector,
                                       Real dx, Real dy, Real h,
                                       Real &newx, Real &newy,
                                       int currentIndex) {
  newx = Xvector[currentIndex] + dx * h;
  newy = Yvector[currentIndex] + dy * h;
}

float laplacian2DGrid::solveLaplace() {

  Real fieldEnergy = 0;
  Real initialValue, newValue;
  for (int v0=1; v0 < this->sizes[0]-1; v0++) {
    for (int v1=1; v1 < this->sizes[1]-1; v1++) {
      initialValue = this->fixedGrid->getVoxel(v0,v1,0);
      if (initialValue > (int)(this->innerValue) &&
          initialValue < (int)(this->outerValue) ) {
        newValue = (this->volume->getVoxel(v0+1, v1, 0) +
                    this->volume->getVoxel(v0, v1+1, 0) + 
                    this->volume->getVoxel(v0-1, v1, 0) + 
                    this->volume->getVoxel(v0, v1-1, 0) ) / 4;
        this->volume->setVoxel(newValue, v0, v1, 0);
        fieldEnergy += sqrt(newValue);
      }
    }
  }
  return fieldEnergy;
}

void laplacian2DGrid::relaxEquation(float convergenceCriteria,
                                    int maxIterations) {

  float convergenceResult1 = 10000; //arbitrary high number
  float convergenceResult2 = 1;     //arbitrary number as well
  float convergence;   //also arbitrary

  int currentIteration = 1;
  
  for (;;) { // loop until break
    convergenceResult1 = this->solveLaplace();

    if (currentIteration > 1) {
      convergence = (convergenceResult2 - convergenceResult1)
	/ convergenceResult1;

      if (this->verbosity >= 1) {
	cout << "Converging: iteration # " << currentIteration 
	     << " with convergence result of " << convergence << endl;
      }
      if ( convergence < convergenceCriteria )
	break;
      else if ( currentIteration >= maxIterations )
        break;
    }
    convergenceResult2 = convergenceResult1;
    currentIteration++;
  }
}

void laplacian2DGrid::createGradients() {
  for (int x=1; x < this->sizes[0]-1; x++) {
    for (int y=1; y < this->sizes[1]-1; y++) {
      if (this->fixedGrid->getVoxel(x,y,0) > this->innerValue &&
          this->fixedGrid->getVoxel(x,y,0) < this->outerValue ) {
        this->gradientX->setVoxel((this->volume->getVoxel(x+1,y,0) - 
                                   this->volume->getVoxel(x-1,y,0)) / 2,
                                  x,y,0);
        this->gradientY->setVoxel((this->volume->getVoxel(x,y+1,0) - 
                                   this->volume->getVoxel(x,y-1,0)) / 2,
                                  x,y,0);
      }
      else {
        this->gradientX->setVoxel(0, x,y,0);
        this->gradientY->setVoxel(0, x,y,0);
      }
    }
  }
}

void laplacian2DGrid::normaliseGradients() {
  Real nx, ny, dx, dy;
  for (int x=1; x < this->sizes[0]-1; x++) {
    for (int y=1; y < this->sizes[1]-1; y++) {
      if (this->fixedGrid->getVoxel(x,y,0) != this->innerValue &&
          this->fixedGrid->getVoxel(x,y,0) != this->outerValue) {
        dx = this->gradientX->getVoxel(x,y,0);
        dy = this->gradientY->getVoxel(x,y,0);

        nx = dx / sqrt( pow(dx,2) + pow(dy,2) );
        ny = dy / sqrt( pow(dy,2) + pow(dx,2) );

        this->gradientX->setVoxel(nx, x,y,0);
        this->gradientY->setVoxel(ny, x,y,0);
      }
    }
  }
}

void laplacian2DGrid::createStreamline(Real x0, Real y0, Real h,
                                       vector<Real> &Xvector,
                                       vector<Real> &Yvector,
                                       interpolation evalType) {
  int i=0;

  Real newx, newy;
  Real dx, dy;

  Xvector.push_back (x0);
  Yvector.push_back (y0); 

  // create vector iterators
  vector<Real>::iterator xIt = Xvector.begin();
  vector<Real>::iterator yIt = Yvector.begin();

  Real evaluation = this->evaluate(Xvector[i],Yvector[i],evalType);

  // move towards outside surface first
  while (evaluation < this->outerValue) {

    this->getDerivatives(Xvector[i], Yvector[i], dx, dy);
    (*this.*integrationStep)(Xvector, Yvector, dx, dy, h,
                             newx, newy, i);
    Xvector.push_back (newx);
    Yvector.push_back (newy);

    i++;
    // omitting test for NaN
    evaluation = this->evaluate(Xvector[i], Yvector[i], evalType);
    // test for runaway resource
    if (Xvector.size() > 50)
      evaluation = this->outerValue;
  }

  // move towards inner surface
  i = 0;

  while (evaluation > this->innerValue) {
    // set the iterators to the beginning of the vectors
    xIt = Xvector.begin();
    yIt = Yvector.begin();

    // always look at the beginning of vector
    this->getDerivatives(Xvector[0], Yvector[0],dx, dy);

    // insert at the beginning, using inverse of h
    (*this.*integrationStep)(Xvector, Yvector, dx, dy, (h * -1),
                             newx, newy, 0);

    Xvector.insert(Xvector.begin(), newx);
    Yvector.insert(Yvector.begin(), newy);

    // omit NaN check for now
    evaluation = this->evaluate(Xvector[0], Yvector[0], evalType);
    i++;
  }
}

Real laplacian2DGrid::streamLength(vector<Real> &Xvector,
                                   vector<Real> &Yvector) {
  
  Real x, y, length;
  length = 0;

  for (int i=0; i < Xvector.size()-1; i++) {
    x = pow(Xvector[i] - Xvector[i+1], 2);
    y = pow(Yvector[i] - Yvector[i+1], 2);
    length += sqrt(x+y);
  }

  // normalise the unit to be millimetres
  Real separations[MAX_DIMENSIONS];
  get_volume_separations (this->fixedGrid->getVolume(), separations );
  // the step size ought to be uniform
  length *= separations[0];
  return length;
}

void laplacian2DGrid::computeAllThickness( Real h, interpolation evalType ) {
  
  
  initialize_progress_report(&this->progressReport, FALSE, this->sizes[0]-1,
                             "Computing thickness streamlines");

  for (int v1=1; v1 < this->sizes[0]-1; v1++) {
    for (int v2=1; v2 < this->sizes[1]-1; v2++) {
      if (this->fixedGrid->getVoxel(v1, v2, 0) > this->innerValue &&
          this->fixedGrid->getVoxel(v1, v2, 0) < this->outerValue) {
        vector<Real> xv, yv;
        this->createStreamline(v1, v2, h, xv, yv, evalType);
        Real length = this->streamLength(xv, yv);
        this->volume->setVoxel(length, v1, v2, 0);
      }
      else {
        this->volume->setVoxel(0, v1, v2, 0);
      }
    }
    update_progress_report( &this->progressReport, v1 );
  }
  terminate_progress_report( &this->progressReport );
}

void laplacian2DGrid::output( char *filename, bool isTextFile ) {
  if ( isTextFile == false ) {
    this->volume->output(filename);
  }
  else {
    cerr << "Text file not supported in 2D mode" << endl;
  }
}
