#ifndef _JPL_POINT_MAP_
#define _JPL_POINT_MAP_

extern "C" {
#include "volume_io.h"
#include "bicpl.h"
}

#include <map>
#include "jplPoint.h"

typedef map<jplPoint, short>::const_iterator pointIterator;

class jplPointMap {
 private:
  map<jplPoint, short> pointMap;
  map<jplPoint, short>::iterator found;
 public:
  jplPointMap();// default constructor
  void addPoint(short coord1, short coord2, short coord3, short value);
  void addPoint(jplPoint point, short value);
  short findPoint(short coord1, short coord2, short coord3);
  short findPoint(jplPoint point);
  int copyToVolume(Volume volume);
  virtual ~jplPointMap();
};

#endif // _JPL_POINT_MAP_
