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

  // and construct the gradient volumes - as copies of volume
  this->gradientX = new mniVolume(this->volume);
  this->gradientY = new mniVolume(this->volume);
  this->gradientZ = new mniVolume(this->volume);

  // set their real ranges - safe here since I don't care about
  // the input volume
  Real lower = -10000;
  Real upper = 10000;
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
  Real value;
  Real topNumerator, bottomNum1, bottomNum2, denom1;

  // create the gradients using central difference formula
  for (int z=1; z < this->sizes[0]-1; z++) {
    for (int x=1; x < this->sizes[1]-1; x++) {
      for (int y=1; y < this->sizes[2]-1; y++) {
        this->gradientZ->setVoxel((this->volume->getVoxel(z+1,x,y) - 
                                   this->volume->getVoxel(z-1,x,y)) / 2,
                                  z,x,y);
        this->gradientX->setVoxel((this->volume->getVoxel(z,x+1,y) - 
                                   this->volume->getVoxel(z,x-1,y)) / 2,
                                  z,x,y);
        this->gradientY->setVoxel((this->volume->getVoxel(z,x,y+1) - 
                                   this->volume->getVoxel(z,x,y-1)) / 2,
                                  z,x,y);
        /*
        if (this->gradientX->getVoxel(z,x,y) < 0) {
          cout << (this->volume->getVoxel(z,x+1,y) - 
                   this->volume->getVoxel(z,x-1,y)) / 2 
              << " "
              << this->gradientX->getVoxel(z,x,y) << endl;
        }
        */
        // normalise the X gradient at that position

        topNumerator = this->gradientX->getVoxel(z,x,y);
        bottomNum1 = this->gradientX->getVoxel(z,x,y);
        bottomNum2 = this->gradientY->getVoxel(z,x,y);
        denom1 = this->gradientZ->getVoxel(z,x,y);

        bottomNum1 = pow(bottomNum1,2);
        bottomNum2 = pow(bottomNum2,2);
        denom1 = pow(denom1, 2);

        value = topNumerator / ((bottomNum1 + bottomNum2) / denom1);

        /*
          value = this->gradientX->getVoxel(z,x,y) / 
          sqrt(((pow(this->gradientX->getVoxel(z,x,y),2)) + 
          pow(this->gradientY->getVoxel(z,x,y),2)) / 
          (pow(this->gradientZ->getVoxel(z,x,y),2)));
        */
        this->gradientX->setVoxel(value, z,x,y);
      }
    }
  }
  this->gradientX->output("gradientX.mnc");
}


void laplacianGrid::output(char *filename) {
  this->volume->output(filename);
}
    

