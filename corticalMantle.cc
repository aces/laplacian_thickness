/*
 * Methods of the corticalMantle Class
 *
 */

#include "corticalMantle.h"

corticalMantle::corticalMantle( STRING outerSurfaceFile,
                                STRING innerSurfaceFile,
                                STRING outputFile,
                                STRING clsFile ) {

  this->outputFile = outputFile;
  this->clsFile = clsFile;
  this->innerSurface = innerSurfaceFile;
  this->outerSurface = outerSurfaceFile;

  // NOTE: this should be made more flexible!
  this->dimNames[0] = MIzspace;
  this->dimNames[1] = MIxspace;
  this->dimNames[2] = MIyspace;

  // get the header information from the cls file - used to get
  // the dimensions for the scanning of the objects
  if (input_volume_header_only(this->clsFile, 3, this->dimNames,
                               &this->mantle, NULL) != OK) {
    cerr << "ERROR: could not open " << this->clsFile << endl;
    exit(1);
  }
  this->mantle = create_label_volume(this->mantle, NC_BYTE);
  get_volume_sizes(this->mantle, this->sizes);

}

void corticalMantle::scanObjectsToVolume(Real maxDistance=1.0,
                                         int innerValue=1,
                                         int outerValue=2) {

  Volume          inner, outer;
  File_formats    format;
  int             obj, n_objects;
  object_struct   **objectsInner, **objectsOuter;
  int            value;

  // set the class variables
  this->greyValue = outerValue;
  this->whiteValue = innerValue;
  this->overlapValue = outerValue + innerValue;

  // create the two label volumes
  inner = create_label_volume(this->mantle, NC_BYTE);
  outer = create_label_volume(this->mantle, NC_BYTE);
  set_all_volume_label_data(inner, 0);
  set_all_volume_label_data(outer, 0);

  // open the inner surface
  if (input_graphics_file(this->innerSurface, &format, &n_objects,
                          &objectsInner) != OK) {
    cerr << "ERROR: could not open " << this->innerSurface << endl;
    exit(1);
  }

  // scan it to the volume
  for_less(obj, 0, n_objects) {
    scan_object_to_volume(objectsInner[obj], this->mantle, inner, 
                          innerValue, maxDistance);
  }

  // open the outer surface
  if (input_graphics_file(this->outerSurface, &format, &n_objects,
                          &objectsOuter) != OK) {
    cerr << "ERROR: could not open " << this->outerSurface << endl;
    exit(1);
  }

  // scan it to the volume
  for_less(obj, 0, n_objects) {
    scan_object_to_volume(objectsOuter[obj], this->mantle, outer, 
                          outerValue, maxDistance);
  }

  // now add the two together
  for (int z=0; z < this->sizes[0]; z++) {
    for (int x=0; x < this->sizes[1]; x++) {
      for (int y=0; y < this->sizes[2]; y++) {
        value = get_volume_label_data_5d(inner, z, x, y, 0, 0) + 
          get_volume_label_data_5d(outer, z, x, y, 0, 0);
        set_volume_label_data_5d(this->mantle, z, x, y, 0, 0, value);
      }
    }
  }

}

int corticalMantle::neighbourFill( int fillValue ) {
  int numValuesChanged = 0;
  int indices[3], tmpIndices[3];

  // loop over the volume, checking for neighbours
  for (indices[0]=1; indices[0] < this->sizes[0]-1; indices[0]++) {
    for (indices[1]=1; indices[1] < this->sizes[1]-1; indices[1]++) {
      for (indices[2]=1; indices[2] < this->sizes[2]-1; indices[2]++) {
        //check if voxel has requisite value
        if (get_volume_label_data(this->mantle, indices) == fillValue ) {
          //check all neighbours
          for (int i=0; i < 3; i++) {
            for (int j=-1; j < 2; j++) {
              tmpIndices = indices;
              tmpIndices[i] += j;
              int value = get_volume_label_data(this->mantle, tmpIndices);
              if (value != this->greyValue &&
                  value != this->whiteValue &&
                  value != this->overlapValue &&
                  value != fillValue) {
                set_volume_label_data(this->mantle, tmpIndices, fillValue);
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

  // fill the outer area first, modifying the mantle volume in place
  // set the initial voxel:
  set_volume_label_data_5d(this->mantle, 5, 5, 5, 0, 0, outerValue);
  
  // recurse through the volume, setting anything which has a neighbour
  // with the outer value but is not labeled as anything else as outer
  // mantle area as well.
  int numValuesChanged = 1;
  while ( numValuesChanged > 0 ) {
    numValuesChanged = this->neighbourFill(outerValue);
    cout << "Changed: " << numValuesChanged <<  endl;
  }

  // now for the inner part
  set_volume_label_data_5d(this->mantle, 63, 108, 105, 0, 0, innerValue);
  numValuesChanged = 1;
  while ( numValuesChanged > 0 ) {
    numValuesChanged = this->neighbourFill(innerValue);
    cout << "Inner: " << numValuesChanged << endl;
  }

  // now fill in the mantle
  for (int v1=0; v1 < this->sizes[0]; v1++) {
    for (int v2=0; v2 < this->sizes[1]; v2++) {
      for (int v3=0; v3 < this->sizes[2]; v3++) {
        int value = get_volume_label_data_5d(this->mantle, v1, v2, v3, 0, 0);
        if (value != innerValue && value != outerValue)
          set_volume_label_data_5d(this->mantle, v1, v2, v3, 0, 0, 
                                   mantleValue);
            
      }
    }
  }    


}

  
            
