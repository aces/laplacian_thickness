#include "jplPointMap.h"

jplPointMap :: jplPointMap() {
  // does nothing - default constructor
}

void jplPointMap :: addPoint(short coord1, short coord2, 
			     short coord3, short value) {
  pointMap[jplPoint (coord1, coord2, coord3)] = value;
}

void jplPointMap :: addPoint(jplPoint point, short value) {
  pointMap[point] = value;
}

short jplPointMap :: findPoint(short coord1, short coord2, short coord3) {
  found = pointMap.find(jplPoint (coord1, 
				  coord2, 
				  coord3));
  return found->second;
}

short jplPointMap :: findPoint(jplPoint point) {
  found = pointMap.find(point);
  return found->second;
}

int jplPointMap :: copyToVolume(Volume volume) {
  int numPointsCopied = 0;
  for (pointIterator m = pointMap.begin(); m != pointMap.end(); ++m) {
    numPointsCopied++;
    // needs check to make sure that none of the points fall outside
    // allowed range for the volume
    jplPoint tmpPoint = m->first;
    set_volume_label_data_5d(volume, tmpPoint.v1, tmpPoint.v2,
			  tmpPoint.v3, 0, 0, m->second);
  }
  return numPointsCopied;
}

jplPointMap :: ~jplPointMap() {
  pointMap.clear();
}

