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

class corticalMantle : bicLabelVolume {
private:
  //! Hold the filenames of the two surfaces
  STRING outerSurface, innerSurface;
  //! Holds the file after which to model the output
  STRING clsFile;
  //! The values used in scanning the objects
  int greyValue, whiteValue, overlapValue;
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
  //! Create laplacian grid
  void initialiseLaplacianGrid( int outerValue,
                                int innerValue,
                                int solvableValue );
  //! Scan objects to volume
  /*!
    Creates an initial volume containing representations of the two 
    surfaces by scanning the object files to a volume.

    Code adapted from David MacDonald's scan_object_to_volume.c
  */
  void scanObjectsToVolume(Real maxDistance=1.0,
                           int innerValue=1,
                           int outerValue=2);
};

#endif // __CORTICALMANTLE__
