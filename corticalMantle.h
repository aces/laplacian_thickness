/*
 * corticalMantle.h - The corticalMantle class
 *
 */

#ifndef __CORTICALMANTLE__
#define __CORTICALMANTLE__

// BIC includes
extern "C" {
#include "volume_io.h"
#include "bicpl.h"
}
#include "mniLabelVolume.h"

// Regular includes
#include <iostream.h>
extern "C" {
#include <unistd.h>
}

//! Class to volumetrically deal with grey and white surface
/*!
  A class desigend to deal volumetrically with the inner and outer
  cortical surface. Specifically designed to initialise a grid for
  solving using the Laplacian thickness metric, though it ought to
  be extendable to more general uses as well
*/

class corticalMantle : public mniLabelVolume {
private:
  //! Hold the filenames of the two surfaces
  STRING outerSurface, innerSurface;
  //! Holds the file after which to model the output
  STRING clsFile;
  //! The values used in scanning the objects
  int greyValue, whiteValue, overlapValue;
  //! The verbosity Level
  int verbosityLevel;
  //! Area fill based on neighbourhood information
  int neighbourFill( int fillValue );
public:
  //! Constructor
  /*!
    Class Constructor
    \param outerSurfaceFile The obj file for the outer surface
    \param innerSurfaceFile The obj file for the inner surface
    \param outputFile The mnc file to place output into
  */
  corticalMantle( STRING outerSurfaceFile,
                  STRING innerSurfaceFile,
                  STRING outputFile,
                  STRING clsFile);

  virtual ~corticalMantle();
  //! Create laplacian grid
  void initialiseLaplacianGrid( int outerValue,
                                int innerValue,
                                int solvableValue );
  //! Scan objects to volume
  /*!
    Creates an initial volume containing representations of the two 
    surfaces by scanning the object files to a volume.

    Code adapted from David MacDonald's scan_object_to_volume.c

    \bug The maxDistance variable is part of the function calls in bicpl,
    but I don't think that it actually is ever used anywhere.
  */
  void scanObjectsToVolume(Real maxDistance=1.0,
                           int innerValue=1,
                           int outerValue=2);
  //! Sets all voxels to the same value
  void setAllVoxels(int value) { 
    set_all_volume_label_data(this->volume, value); };
  //! Sets verbosity
  /*!
    Sets the amount of info to be spewed out during processing

    \param verbosityLevel The higher the level, the more will be spewed out.
    - 0 tells it to be quiet. This is the default
    - 1 prints out a general progress report
    - >1 prints out a slew of mostly useless info
  */
  void setVerbosity(int verbosityLevel) {
    this->verbosityLevel = verbosityLevel; };
  

};

#endif // __CORTICALMANTLE__
