#include "laplacianGrid.h"


// void laplacianGrid::*integrationStep(vector<Real> &Xvector,
// 				     vector<Real> &Yvector,
// 				     vector<Real> &Zvector,
// 				     Real dx, Real dy, Real dz, Real h,
// 				     vector<Real>::iterator XinsertIt,
// 				     vector<Real>::iterator YinsertIt,
// 				     vector<Real>::iterator ZinsertIt,
// 				     int currentStep) {}

inline void laplacianGrid::getDerivatives( Real x, Real y, Real z,
					   Real &dx, Real &dy, Real &dz ) {

  dx = this->gradientX->getInterpolatedVoxel(x,y,z,2);
  dy = this->gradientY->getInterpolatedVoxel(x,y,z,2);
  dz = this->gradientZ->getInterpolatedVoxel(x,y,z,2);
}

// take an integration step using Euler's method
inline void laplacianGrid::eulerStep(vector<Real> &Xvector, 
				     vector<Real> &Yvector, 
				     vector<Real> &Zvector,
				     Real dx, Real dy, Real dz, Real h, 
				     vector<Real>::iterator XinsertIt,
				     vector<Real>::iterator YinsertIt,
				     vector<Real>::iterator ZinsertIt,
				     int currentIndex) {

  Xvector.insert( XinsertIt, Xvector[currentIndex] + dx * h );
  Yvector.insert( YinsertIt, Yvector[currentIndex] + dy * h );
  Zvector.insert( ZinsertIt, Zvector[currentIndex] + dz * h );

}

inline void laplacianGrid::secondOrderRungeKuttaStep(
					  vector<Real> &Xvector, 
					  vector<Real> &Yvector, 
					  vector<Real> &Zvector,
					  Real dx, Real dy, Real dz, Real h, 
					  vector<Real>::iterator XinsertIt,
					  vector<Real>::iterator YinsertIt,
					  vector<Real>::iterator ZinsertIt,
					  int currentIndex) {

  Real hh = h * 0.5;
  Real kx, ky, kz;
  Real ddx, ddy, ddz;

  kx = Xvector[currentIndex] + dx * hh;
  ky = Yvector[currentIndex] + dy * hh;
  kz = Zvector[currentIndex] + dz * hh;
						
  this->getDerivatives(kx, ky, kz, ddx, ddy, ddz);

  Xvector.insert( XinsertIt, Xvector[currentIndex] + ddx * h );
  Yvector.insert( YinsertIt, Yvector[currentIndex] + ddy * h );
  Zvector.insert( ZinsertIt, Zvector[currentIndex] + ddz * h );
}
  

// take an integration step using a fourth order Runge-Kutta Model
inline void laplacianGrid::fourthOrderRungeKuttaStep(vector<Real> &Xvector, 
					  vector<Real> &Yvector, 
					  vector<Real> &Zvector,
					  Real dx, Real dy, Real dz, Real h, 
					  vector<Real>::iterator XinsertIt,
					  vector<Real>::iterator YinsertIt,
					  vector<Real>::iterator ZinsertIt,
					  int currentIndex) {

  Real ddx[3], ddy[3], ddz[3];
  Real kx, ky, kz;

  Real hh = h*0.5;

  // take first step
  kx = Xvector[currentIndex] + dx * hh;
  ky = Yvector[currentIndex] + dy * hh;
  kz = Zvector[currentIndex] + dz * hh;

  // second step
  this->getDerivatives(kx, ky, kz, ddx[0], ddy[0], ddz[0]);

  kx = Xvector[currentIndex] + ddx[0] * hh;
  ky = Yvector[currentIndex] + ddy[0] * hh;
  kz = Zvector[currentIndex] + ddz[0] * hh;
  
  // third step
  this->getDerivatives(kx, ky, kz, ddx[1], ddy[1], ddz[1]);
  
  // fourth step
  kx = Xvector[currentIndex] + ddx[1] * h;
  ky = Yvector[currentIndex] + ddy[1] * h;
  kz = Zvector[currentIndex] + ddz[1] * h;

  this->getDerivatives(kx, ky, kz, ddx[2], ddy[2], ddz[2]);
  
  Xvector.insert( XinsertIt, Xvector[currentIndex] + 
		  (dx / 6) + (ddx[0] / 3) + (ddx[1]) / 3 + (ddx[2] / 6) );
  Yvector.insert( YinsertIt, Yvector[currentIndex] +
		  (dy / 6) + (ddy[0] / 3) + (ddy[1]) / 3 + (ddy[2] / 6) );
  Zvector.insert( ZinsertIt, Zvector[currentIndex] + 
		  (dz / 6) + (ddz[0] / 3) + (ddz[1]) / 3 + (ddz[2] / 6) );

}


// constructor from file
laplacianGrid::laplacianGrid(char* mantleFile, 
                             int innerValue, 
                             int outerValue,
			     integrator type) {

  this->innerValue = innerValue;
  this->outerValue = outerValue;

  // create the grid volume from file
  this->fixedGrid = new mniVolume(mantleFile);

  // construct the volume from the mantle file, but signed
  this->volume = new mniVolume(mantleFile,
                               outerValue * -1,
                               outerValue,
                               3,
                               XYZdimOrder,
                               NC_SHORT,
                               TRUE,
                               TRUE,
                               NULL);
  
  this->initialiseVolumes(type);
}

// constructor from volume_struct
laplacianGrid::laplacianGrid(Volume mantleVolume,
			     int innerValue,
			     int outerValue,
			     integrator type) {
  this->innerValue = innerValue;
  this->outerValue = outerValue;

  // take the volume_struct pointer as the grid
  this->fixedGrid = new mniVolume(mantleVolume);
  this->volume = new mniVolume(this->fixedGrid, FALSE, NC_SHORT, TRUE);

  this->initialiseVolumes(type);
}
    

void laplacianGrid::initialiseVolumes(integrator type) {

  // and construct the gradient volumes - using volume definition copy
  this->gradientX = new mniVolume(this->volume, FALSE, NC_SHORT, TRUE, -1, 1);
  this->gradientY = new mniVolume(this->volume, FALSE, NC_SHORT, TRUE, -1, 1);
  this->gradientZ = new mniVolume(this->volume, FALSE, NC_SHORT, TRUE, -1, 1);


  // set their real ranges - safe here since I don't care about
  // the input volume
  Real lower = -1;
  Real upper = 1;
  gradientX->setRealRange(lower, upper);
  gradientY->setRealRange(lower, upper);
  gradientZ->setRealRange(lower, upper);

  // initialise the sizes member
  this->sizes = new int[3];
  this->sizes = this->fixedGrid->getSizes();

  // set default verbosity to 0
  this->verbosity = 0;
  
  // set the function pointer
  cout << "Using integration type: " << type << endl;
  if (type == EULER)
    this->integrationStep = &laplacianGrid::eulerStep;
  else if (type == SECOND_ORDER_RK)
    this->integrationStep = &laplacianGrid::secondOrderRungeKuttaStep;
  else if (type == FOURTH_ORDER_RK)
    this->integrationStep = &laplacianGrid::fourthOrderRungeKuttaStep;

}

// one iteration of solving laplace's equation
float laplacianGrid::solveLaplace() {

  Real fieldEnergy = 0;
  Real initialValue, newValue;
  for (int v0=1; v0 < this->sizes[0]-1; v0++) {
    for (int v1=1; v1 < this->sizes[1]-1; v1++) {
      for (int v2=1; v2 < this->sizes[2]-1; v2++) {
        initialValue = this->fixedGrid->getVoxel(v0,v1,v2);
	//	cout << initialValue << endl;
        if (initialValue > (this->innerValue) &&
            initialValue < (this->outerValue) ) {
	  //cout << "indeed" << endl;
          newValue = (this->volume->getVoxel(v0+1, v1, v2) + 
                      this->volume->getVoxel(v0, v1+1, v2) +
                      this->volume->getVoxel(v0, v1, v2+1) +
                      this->volume->getVoxel(v0-1, v1, v2) +
                      this->volume->getVoxel(v0, v1-1, v2) +
                      this->volume->getVoxel(v0, v1, v2-1) ) / 6;
          this->volume->setVoxel(newValue, v0, v1, v2);
          fieldEnergy += sqrt(newValue);
	  //	  cout << fieldEnergy << endl;
        }
      }
    }
  }

  return fieldEnergy; //should add code to evaluate convergence
}
        
// relax the equation
void laplacianGrid::relaxEquation(float convergenceCriteria,
                                  int maxIterations) {

  float convergenceResult1 = 10000; //arbitrary high number
  float convergenceResult2 = 1;     //arbitrary number as well
  float convergence;   //also arbitrary

  int currentIteration = 1;

  //  while (convergence > convergenceCriteria &&
  //         currentIteration < maxIterations) {
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

// create a gradient volume in each of the cardinal directions
void laplacianGrid::createGradients() {

  // create the gradients using central difference formula,
  // e.g. dx = volume(x+1,y,z) - volume(x-1,y,z,) / 2
  for (int x=1; x < this->sizes[0]-1; x++) {
    for (int y=1; y < this->sizes[1]-1; y++) {
      for (int z=1; z < this->sizes[2]-1; z++) {
        if (this->fixedGrid->getVoxel(x,y,z) > this->innerValue &&
            this->fixedGrid->getVoxel(x,y,z) < this->outerValue ) {
          this->gradientX->setVoxel((this->volume->getVoxel(x+1,y,z) - 
                                     this->volume->getVoxel(x-1,y,z)) / 2,
                                    x,y,z);
          this->gradientY->setVoxel((this->volume->getVoxel(x,y+1,z) - 
                                     this->volume->getVoxel(x,y-1,z)) / 2,
                                    x,y,z);
          this->gradientZ->setVoxel((this->volume->getVoxel(x,y,z+1) - 
                                     this->volume->getVoxel(x,y,z-1)) / 2,
                                    x,y,z);
          
        }
        else {
          this->gradientX->setVoxel(0, x,y,z);
          this->gradientY->setVoxel(0, x,y,z);
          this->gradientZ->setVoxel(0, x,y,z);
        }
      }
    }
  }
}

void laplacianGrid::normaliseGradients() {
	  /* normalise the X gradient at that position like so:
	   * Nx = dx / [dx^2 + dy^2 / dz^2]^0.5
	   */
  Real nx, ny, nz, dx, dy, dz;

  for (int x=1; x < this->sizes[0]-1; x++) {
    for (int y=1; y < this->sizes[1]-1; y++) {
      for (int z=1; z < this->sizes[2]-1; z++) {

        if (this->fixedGrid->getVoxel(x,y,z) != this->innerValue &&
            this->fixedGrid->getVoxel(x,y,z) != this->outerValue) {
  
          dx = this->gradientX->getVoxel(x,y,z);
          dy = this->gradientY->getVoxel(x,y,z);
          dz = this->gradientZ->getVoxel(x,y,z);

          nx = dx / sqrt( pow(dx,2) + pow(dy,2) / pow(dz,2) );
          ny = dy / sqrt( pow(dy,2) + pow(dx,2) / pow(dz,2) );
          nz = dz / sqrt( pow(dz,2) + pow(dy,2) / pow(dx,2) );
	  //          cout << "NX: " << nx << endl;
	  //          cout << "NY: " << ny << endl;
	  //          cout << "NZ: " << nz << endl;
	  
          this->gradientX->setVoxel(nx, x,y,z);
          this->gradientY->setVoxel(ny, x,y,z);
          this->gradientZ->setVoxel(nz, x,y,z);
        }
      }
    }
  }

  //  cout << "Test: " << this->gradientX->getVoxel(123,96,108) << endl;
  //  this->gradientX->output("gradientX.mnc");
  //  this->gradientY->output("gradientY.mnc");
  //  this->gradientZ->output("gradientZ.mnc");
}

    

// use Eulers method, first towards one then towards the other surface
void laplacianGrid::createStreamline(int x0, int y0, int z0, Real h,
                                     vector<Real> &Xvector, 
                                     vector<Real> &Yvector,
                                     vector<Real> &Zvector) {

  int i = 0;

  //  Xvector[i] = (x0); 
  Xvector.push_back (x0);
  Yvector.push_back (y0); 
  Zvector.push_back (z0); 

  // create vector iterators
  vector<Real>::iterator xIt = Xvector.begin();
  vector<Real>::iterator yIt = Yvector.begin();
  vector<Real>::iterator zIt = Zvector.begin();

  Real evaluation = this->fixedGrid->getInterpolatedVoxel(Xvector[i],
							  Yvector[i],
							  Zvector[i]); 


  if (this->verbosity >= 5) {
    cout << "Start of createStreamline function." << endl;
    cout << "First evaluation: " << evaluation << endl;
    cout << "Initial vectors (xyz): " << Xvector[i] << " "
         << Yvector[i] << " " << Zvector[i] << endl;
    cout << "Using H of: " << h << endl << endl;
  }

  // move towards outside surface first
  while (evaluation < this->outerValue) {

    //    cout << Xvector[i] << " " << Yvector[1] << " " << Zvector[i] << endl;

    Real Xvalue = this->gradientX->getInterpolatedVoxel(Xvector[i],
                                                        Yvector[i],
                                                        Zvector[i],
                                                        0);
    Real Yvalue = this->gradientY->getInterpolatedVoxel(Xvector[i],
                                                        Yvector[i],
                                                        Zvector[i],
                                                        0);
    Real Zvalue = this->gradientZ->getInterpolatedVoxel(Xvector[i],
                                                        Yvector[i],
                                                        Zvector[i],
                                                        0);

    //    Xvector.push_back (Xvector[i] + Xvalue * h);
    //    Yvector.push_back (Yvector[i] + Yvalue * h);
    //    Zvector.push_back (Zvector[i] + Zvalue * h);

    (*this.*integrationStep)(Xvector, Yvector, Zvector, 
			  Xvalue, Yvalue, Zvalue, h,
			  Xvector.end(), Yvector.end(), Zvector.end(), i);

    //eulerStep(Xvector, Yvector, Zvector, Xvalue, Yvalue, Zvalue, h,
    //      Xvector.end(), Yvector.end(), Zvector.end(), 
    //      xIt, yIt, zIt);
	      

     if (this->verbosity >= 5) {
       cout << "Interpolated values: " << Xvalue << " " 
            << Yvalue << " " << Zvalue << endl;
       cout << "Values pushed back: " << Xvector[i] + Xvalue * h  << " "
            << Yvector[i] + Yvalue * h << " "
            << Zvector[i] + Zvalue * h << endl;
       cout << "I = " << i << endl << endl;
     }

     //    cout << i << endl;
    i++;
    //    xIt++; yIt++; zIt++;
    // test for Not A Number
    if (Zvector[i] != Zvector[i] || 
        Xvector[i] != Xvector[i] || 
        Yvector[i] != Yvector[i]) {
      if (this->verbosity >= 2) {
        cerr << "WARNING: NaN at xyz: " << x0 << " " << y0 << " " 
             << z0 << endl;
        cerr << "with values: " << Xvector[i] << " " << Yvector[i] << " "
             << Zvector[i] << endl;
      }
      // end search
      return;
    }
    else {
      evaluation = this->fixedGrid->getInterpolatedVoxel((Xvector[i]),
							 (Yvector[i]),
							 (Zvector[i]));
    
      if (this->verbosity >=5 ) {
        cout << "Evaluation function: " << evaluation << endl << endl;
      }
      // test for runaway resource
      if (Xvector.size() > 50)
        evaluation = this->outerValue;
    }
    if (this->verbosity >= 5) {
      cout << "--------------------------" << endl << endl;
    }
  }

  if (this->verbosity >= 5) {
    cout << "Starting descent towards inner surface." << endl;
  }

  // move towards inner surface
  while (evaluation != this->innerValue) {
    // set the iterators to the beginning of the vectors
    xIt = Xvector.begin();
    yIt = Yvector.begin();
    zIt = Zvector.begin();

    // test for NaN
    // always look at the beginning of vector
    Real Xvalue = this->gradientX->getInterpolatedVoxel(Xvector[0],
                                                        Yvector[0],
                                                        Zvector[0],
                                                        2);
    Real Yvalue = this->gradientY->getInterpolatedVoxel(Xvector[0],
                                                        Yvector[0],
                                                        Zvector[0],
                                                        2);
    Real Zvalue = this->gradientZ->getInterpolatedVoxel(Xvector[0],
                                                        Yvector[0],
                                                        Zvector[0],
                                                        2);
    
    // insert at the beginning, using inverse of h
    (*this.*integrationStep)(Xvector, Yvector, Zvector, 
			  Xvalue, Yvalue, Zvalue, (h * -1),
			  Xvector.begin(), Yvector.begin(), Zvector.begin(), 
			  0);

//     Xvector.insert(xIt, Xvector[0] + Xvalue * (h * -1));
//     Yvector.insert(yIt, Yvector[0] + Yvalue * (h * -1));
//     Zvector.insert(zIt, Zvector[0] + Zvalue * (h * -1));

    if (this->verbosity >= 5) {
      cout << "Interpolated values: " << Xvalue << " " 
           << Yvalue << " " << Zvalue << endl;
      cout << "Values pushed back: " << Xvector[0] + Xvalue * h  << " "
           << Yvector[0] + Yvalue * h << " "
           << Zvector[0] + Zvalue * h << endl;
      cout << "I = " << i << endl << endl;
    }

    
    if (Zvector[0] != Xvector[0] || 
        Xvector[0] != Yvector[0] || 
        Yvector[0] != Zvector[0]) {
      if (this->verbosity >=2 ) {
        cerr << "WARNING: NaN at xyz: " << x0 << " " << y0 << " " 
             << z0 << endl;
      }
      return;
    }
    else {

      evaluation = this->fixedGrid->getVoxel((int)rint(Xvector[0]),
                                             (int)rint(Yvector[0]),
                                             (int)rint(Zvector[0]));
      if (Xvector.size() > 50)
        evaluation = this->innerValue;
      
    }
    if (this->verbosity >= 5) {
      cout << "--------------------------" << endl << endl;
    }
  }
  if (this->verbosity >= 5) {
    cout << "End of createStreamline function" << endl;
  }
}
  
Real laplacianGrid::streamLength(vector<Real> &Xvector,
                                 vector<Real> &Yvector,
                                 vector<Real> &Zvector) {

  Real x, y, z, length;
  length = 0;

  for (int i=0; i < Xvector.size()-1; i++) {
    x = pow(Xvector[i] - Xvector[i+1], 2);
    y = pow(Yvector[i] - Yvector[i+1], 2);
    z = pow(Zvector[i] - Zvector[i+1], 2);
    length += sqrt(x + y + z);
  }
  
  // normalise unit to be millimetres
  Real separations[MAX_DIMENSIONS];
  get_volume_separations( this->fixedGrid->getVolume(), separations );
  // the step size ought to be uniform. If not we're screwed anyway.
  length *= separations[0];
  return length;
}

void laplacianGrid::computeAllThickness( Real h, char *objFile ) {
	File_formats 	format;
	int 		point, nPoints, nObjects;
	Point 		*points;
	object_struct 	**objects;
	Real            *voxel, length;
		
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
	
	initialize_progress_report( &this->progressReport, FALSE, nPoints,
	                            "Computing thickness streamlines" );
 	for (int i=0; i < nPoints; i++) {
 	  voxel = this->fixedGrid->convertWorldToVoxel(
 	                            RPoint_x(points[i]),
 	                            RPoint_y(points[i]),
 	                            RPoint_z(points[i]) );
 	
  	vector<Real> xv, yv, zv;
//  	cout << voxel[0] << " " << voxel[1] << " " << voxel[2] << endl;
  	if (this->fixedGrid->getVoxel(int(rint(voxel[0])), 
				      int(rint(voxel[1])), 
				      int(rint(voxel[2]))) 
				      < 6000 &&
				      this->fixedGrid->getVoxel
				      (int(rint(voxel[0])), 
				       int(rint(voxel[1])), 
				       int(rint(voxel[2]))) 
				      > 4000) {
  	  this->createStreamline( (int)rint(voxel[0]),
	                            (int)rint(voxel[1]),
	                            (int)rint(voxel[2]),
	                            h, xv, yv, zv );
    	length = this->streamLength( xv, yv, zv );
//    	cout << "Actually computing streamlines." << endl;
   } else {
      length = 0;
      if ( this->verbosity > 1 ) {
        cout << "Could not evaluate at vertex " << i << endl;
      }
   }
	  this->thicknessPerVertex[i] = length;
	  update_progress_report( &this->progressReport, i );
	}
}
	

void laplacianGrid::computeAllThickness(Real h) {

  this->volume->setRealRange(0, 100);

  initialize_progress_report(&this->progressReport, FALSE, this->sizes[0]-1,
                             "Computing thickness streamlines");

  for (int v1=1; v1 < this->sizes[0]-1; v1++) {
    for (int v2=1; v2 < this->sizes[1]-1; v2++) {
      for (int v3=1; v3 < this->sizes[2]-1; v3++) {
        if (this->fixedGrid->getVoxel(v1, v2, v3) > this->innerValue && 
            this->fixedGrid->getVoxel(v1, v2, v3) < this->outerValue) {
	  //	  cout << v1 << " " << v2 << " " << v3 << endl;
          vector<Real> xv, yv, zv;
          this->createStreamline(v1, v2, v3, h, xv, yv, zv);
          Real length = this->streamLength(xv, yv, zv);
          this->volume->setVoxel(length ,v1, v2, v3);
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



    


