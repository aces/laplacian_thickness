#include "laplacianGrid.h"

// constructor from file
laplacianGrid::laplacianGrid(char* mantleFile, 
                             int innerValue, 
                             int outerValue) {

  this->innerValue = innerValue;
  this->outerValue = outerValue;

  // create the grid volume from file
  this->fixedGrid = new mniLabelVolume(mantleFile);
                                
  
  // construct the volume from the mantle file, but signed
  this->volume = new mniVolume(mantleFile,
                               outerValue * -1,
                               outerValue,
                               3,
                               ZXYdimOrder,
                               NC_SHORT,
                               TRUE,
                               TRUE,
                               NULL);

  // and construct the gradient volumes - using volume definition copy
  this->gradientX = new mniVolume(this->volume, TRUE, NC_FLOAT, TRUE, -1, 1);
  this->gradientY = new mniVolume(this->volume, TRUE, NC_FLOAT, TRUE, -1, 1);
  this->gradientZ = new mniVolume(this->volume, TRUE, NC_FLOAT, TRUE, -1, 1);


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
}

// one iteration of solving laplace's equation
float laplacianGrid::solveLaplace() {

  Real initialValue, newValue;
  for (int v0=1; v0 < this->sizes[0]-1; v0++) {
    for (int v1=1; v1 < this->sizes[1]-1; v1++) {
      for (int v2=1; v2 < this->sizes[2]-1; v2++) {
        initialValue = this->fixedGrid->getVoxel(v0,v1,v2);
        if (initialValue == this->innerValue ||
            initialValue == this->outerValue) {
          // do nothing
        }
        else {
          newValue = (this->volume->getVoxel(v0+1, v1, v2) + 
                      this->volume->getVoxel(v0, v1+1, v2) +
                      this->volume->getVoxel(v0, v1, v2+1) +
                      this->volume->getVoxel(v0-1, v1, v2) +
                      this->volume->getVoxel(v0, v1-1, v2) +
                      this->volume->getVoxel(v0, v1, v2-1) ) / 6;
          this->volume->setVoxel(newValue, v0, v1, v2);
        }
      }
    }
  }

  return 0.0; //should add code to evaluate convergence
}
        
// relax the equation
void laplacianGrid::relaxEquation(float convergenceCriteria,
                                  int maxIterations) {

  float convergenceResult = 10000; //arbitrary high number
  int currentIteration = 1;

  while (convergenceResult > convergenceCriteria &&
         currentIteration < maxIterations) {
    convergenceResult = this->solveLaplace();
    currentIteration++;
  }
}

// create a gradient volume in each of the cardinal directions
void laplacianGrid::createGradients() {

  // create the gradients using central difference formula,
  // e.g. dx = volume(x+1,y,z) - volume(x-1,y,z,) / 2
  for (int z=1; z < this->sizes[0]-1; z++) {
    for (int x=1; x < this->sizes[1]-1; x++) {
      for (int y=1; y < this->sizes[2]-1; y++) {
        if (this->fixedGrid->getVoxel(z,x,y) != this->innerValue &&
            this->fixedGrid->getVoxel(z,x,y) != this->outerValue) {
          this->gradientZ->setVoxel((this->volume->getVoxel(z+1,x,y) - 
                                     this->volume->getVoxel(z-1,x,y)) / 2,
                                    z,x,y);
          this->gradientX->setVoxel((this->volume->getVoxel(z,x+1,y) - 
                                     this->volume->getVoxel(z,x-1,y)) / 2,
                                    z,x,y);
          this->gradientY->setVoxel((this->volume->getVoxel(z,x,y+1) - 
                                     this->volume->getVoxel(z,x,y-1)) / 2,
                                    z,x,y);
          
        }
        else {
          this->gradientX->setVoxel(0, z,x,y);
          this->gradientY->setVoxel(0, z,x,y);
          this->gradientZ->setVoxel(0, z,x,y);
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

  for (int z=1; z < this->sizes[0]-1; z++) {
    for (int x=1; x < this->sizes[1]-1; x++) {
      for (int y=1; y < this->sizes[2]-1; y++) {

        if (this->fixedGrid->getVoxel(z,x,y) != this->innerValue &&
            this->fixedGrid->getVoxel(z,x,y) != this->outerValue) {
  
          dx = this->gradientX->getVoxel(z,x,y);
          dy = this->gradientY->getVoxel(z,x,y);
          dz = this->gradientZ->getVoxel(z,x,y);

          nx = dx / sqrt( pow(dx,2) + pow(dy,2) / pow(dz,2) );
          ny = dy / sqrt( pow(dy,2) + pow(dx,2) / pow(dz,2) );
          nz = dz / sqrt( pow(dz,2) + pow(dx,2) / pow(dy,2) );
	  //          cout << "NX: " << nx << endl;
	  //          cout << "NY: " << ny << endl;
	  //          cout << "NZ: " << nz << endl;
	  
          this->gradientX->setVoxel(nx, z,x,y);
          this->gradientY->setVoxel(ny, z,x,y);
          this->gradientZ->setVoxel(nz, z,x,y);
        }
      }
    }
  }

  cout << "Test: " << this->gradientX->getVoxel(123,96,108) << endl;
  this->gradientX->output("gradientX.mnc");
  this->gradientY->output("gradientY.mnc");
  this->gradientZ->output("gradientZ.mnc");
  cout << "Test: " << this->gradientX->getVoxel(123,96,108) << endl;
}

// use Eulers method, first towards one then towards the other surface
void laplacianGrid::createStreamline(int x0, int y0, int z0, int h,
                                     vector<Real> &Xvector, 
                                     vector<Real> &Yvector,
                                     vector<Real> &Zvector) {

  int i = 0;
  // initialise a new vector
  Xvector[i] = x0; 
  Yvector[i] = y0; 
  Zvector[i] = z0; 


    cout << "I: " << i << " XV: " << Xvector[i] 
	 << " YV: " << Yvector[i] 
	 << " ZV: " << Zvector[i] << endl;

cout << "Test: " << this->gradientX->getVoxel(123,96,108) << endl;
  // move towards outside surface first
  while (this->fixedGrid->getVoxel((int)Zvector[i],
				   (int)Xvector[i],
				   (int)Yvector[i]) 
         != this->outerValue) {
    Real Xvalue = this->gradientX->getInterpolatedVoxel(Zvector[i],Xvector[i],Yvector[i]);
    Real Yvalue = this->gradientY->getInterpolatedVoxel(Zvector[i],Xvector[i],Yvector[i]);
    Real Zvalue = this->gradientZ->getInterpolatedVoxel(Zvector[i],Xvector[i],Yvector[i]);
    cout << Xvalue << " " << Yvalue << " " << Zvalue << endl;
    cout << Xvalue << " " << this->gradientX->getVoxel(123,108,96) << endl;

    Xvector.push_back(Xvector[i] + Xvalue * h);
    Yvector.push_back(Yvector[i] + Yvalue * h);
    Zvector.push_back(Zvector[i] + Zvalue * h);


    cout << "I: " << i << " XV: " << Xvector[i] 
	 << " YV: " << Yvector[i] 
	 << " ZV: " << Zvector[i] 
	 << " Grid " << this->fixedGrid->getVoxel((int)Zvector[i],
						(int)Xvector[i],
						(int)Yvector[i])
	 << " test " << this->gradientX->getVoxel(123,108,96)
	 <<endl;
    i++;
  }
}
  
  
				     

void laplacianGrid::output(char *filename) {
  this->volume->output(filename);
}
    

