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

  
            
