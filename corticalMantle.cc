/*
 * Methods of the corticalMantle Class
 *
 */

#include "corticalMantle.h"

corticalMantle::corticalMantle( STRING outerSurfaceFile,
                                STRING innerSurfaceFile,
                                STRING outputFile,
                                STRING clsFile ) 
  : mniLabelVolume(clsFile, 1, 3, ZXYdimOrder, NC_SHORT, NULL) {

  this->filename = outputFile;
  this->clsFile = clsFile;
  this->innerSurface = innerSurfaceFile;
  this->outerSurface = outerSurfaceFile;

  // set default verbosity to 0, or quiet
  this->verbosityLevel = 0;

}

corticalMantle::~corticalMantle() {
  delete_volume(this->volume);
  delete this->sizes;
}

void corticalMantle::scanObjectsToVolume(Real maxDistance=1.0,
                                         int innerValue=1,
                                         int outerValue=2) {

  mniLabelVolume  *inner, *outer;
  File_formats    format;
  int             obj, n_objects;
  object_struct   **objectsInner, **objectsOuter;
  int            value;

  if (this->verbosityLevel > 0) {
    cout << "Starting the scan of the objects to the volume ..." << endl;
  }

  // set the class variables
  this->greyValue = outerValue;
  this->whiteValue = innerValue;
  this->overlapValue = outerValue + innerValue;

  // create the two label volumes
  inner = new mniLabelVolume(this);
  outer = new mniLabelVolume(this);
  inner->setAllVoxels(0);
  outer->setAllVoxels(0);

  // open the inner surface
  if (input_graphics_file(this->innerSurface, &format, &n_objects,
                          &objectsInner) != OK) {
    cerr << "ERROR: could not open " << this->innerSurface << endl;
    exit(1);
  }

  // scan it to the volume
  if (this->verbosityLevel > 0) {
    cout << "Scanning " << this->innerSurface << " to volume." << endl;
  }
  for_less(obj, 0, n_objects) {
    scan_object_to_volume(objectsInner[obj], this->volume, inner->getVolume(), 
                          innerValue, maxDistance);
  }

  // open the outer surface
  if (input_graphics_file(this->outerSurface, &format, &n_objects,
                          &objectsOuter) != OK) {
    cerr << "ERROR: could not open " << this->outerSurface << endl;
    exit(1);
  }

  // scan it to the volume
  if (this->verbosityLevel > 0) {
    cout << "Scanning " << this->outerSurface << " to volume." << endl;
  }
  for_less(obj, 0, n_objects) {
    scan_object_to_volume(objectsOuter[obj], this->volume, outer->getVolume(), 
                          outerValue, maxDistance);
  }

  // now add the two together
  if (this->verbosityLevel < 0) {
    cout << "Creating final objects scan volume." << endl;
  }
  for (int z=0; z < this->sizes[0]; z++) {
    for (int x=0; x < this->sizes[1]; x++) {
      for (int y=0; y < this->sizes[2]; y++) {
        this->setVoxel(inner->getVoxel(z,x,y) + outer->getVoxel(z,x,y),
                       z,x,y);
      }
    }
  }
  delete inner;
  delete outer;

}

int corticalMantle::neighbourFill( int fillValue ) {
  int numValuesChanged = 0;
  int indices[3], tmpIndices[3];

  // loop over the volume, checking for neighbours

  for (indices[0]=1; indices[0] < this->sizes[0]-1; indices[0]++) {
    for (indices[1]=1; indices[1] < this->sizes[1]-1; indices[1]++) {
      for (indices[2]=1; indices[2] < this->sizes[2]-1; indices[2]++) {
        //check if voxel has requisite value
        //        int tmp = this->getVoxel(indices);
        //        cout << indices[0] << " " << indices[1] << " " << indices[2] 
        //             << " " << tmp << endl;
        if (this->getVoxel(indices) == fillValue ) {
          //check all neighbours
          for (int i=0; i < 3; i++) {
            for (int j=-1; j < 2; j++) {
              tmpIndices = indices;
              tmpIndices[i] += j;
              int value = this->getVoxel(tmpIndices);
              if (value != this->greyValue &&
                  value != this->whiteValue &&
                  value != this->overlapValue &&
                  value != fillValue) {
                this->setVoxel(fillValue, tmpIndices);
                //cout << "Indices: " << tmpIndices[0] << " " << tmpIndices[1]
                //   << " " << tmpIndices[2] << " " << value << " " 
                //   << fillValue << endl;
                numValuesChanged++;
              }
            }
          }
        }
      }
    }
  }
  return numValuesChanged;
}

void corticalMantle::initialiseLaplacianGrid( int outerValue,
                                              int innerValue,
                                              int mantleValue ) {

  // temporary values to use for filling - to avoid conflicts
  int tmpOuterValue = 5;
  int tmpInnerValue = 6;
  int tmpMantleValue = 7;

  // fill the outer area first, modifying the mantle volume in place
  // set the initial voxel:
  this->setVoxel(tmpOuterValue, 5, 5, 5);
  
  // recurse through the volume, setting anything which has a neighbour
  // with the outer value but is not labeled as anything else as outer
  // mantle area as well.
  if (this->verbosityLevel > 0)
    cout << "Beginning relaxation of area outside of cortex. Will stop once a value of 0 voxels changed has been reached." << endl;

  int numValuesChanged = 1;
  while ( numValuesChanged > 0 ) {
    numValuesChanged = this->neighbourFill(tmpOuterValue);
    if (this->verbosityLevel > 0)
      cout << "Voxels changed: " << numValuesChanged <<  endl;
  }

  // now for the inner part
  this->setVoxel(tmpInnerValue, 63, 108, 105);
  numValuesChanged = 1;
  
  if (this->verbosityLevel > 0)
    cout << "Beginning relaxation of area inside of cortex. Will stop once a value of 0 voxels changed has been reached." << endl;
  while ( numValuesChanged > 0 ) {
    numValuesChanged = this->neighbourFill(tmpInnerValue);
    if (this->verbosityLevel > 0)
      cout << "Voxels changed: " << numValuesChanged << endl;
  }

  // now fill in the mantle
  if (this->verbosityLevel > 0)
    cout << "Creating final grid volume." << endl;
  for (int v1=0; v1 < this->sizes[0]; v1++) {
    for (int v2=0; v2 < this->sizes[1]; v2++) {
      for (int v3=0; v3 < this->sizes[2]; v3++) {
        int value = this->getVoxel(v1, v2, v3);
        if (value == tmpInnerValue)
          this->setVoxel(innerValue, v1, v2, v3);
        if (value == tmpOuterValue)
          this->setVoxel(outerValue, v1, v2, v3);
        if (value != tmpInnerValue && value != tmpOuterValue)
          this->setVoxel(mantleValue, v1, v2, v3);
            
      }
    }
  }    


}
