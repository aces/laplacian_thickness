#include "jplPoint.h"

jplPoint :: jplPoint(short point1, short point2, short point3) {
  v1 = point1;
  v2 = point2;
  v3 = point3;
}

jplPoint :: jplPoint() {
  // nothing here
}

void jplPoint ::  setPoints(short point1, short point2, short point3) {
  v1 = point1;
  v2 = point2;
  v3 = point3;
}

ostream& operator<< (ostream& out, const jplPoint& aPoint) {
  return out << aPoint.v1 << " " << aPoint.v2 << " " << aPoint.v3;
}


bool jplPoint :: operator< (const jplPoint& aPoint) const {
  // sorting rule: v1 < v2 < v3 in terms of priorities 
  if (aPoint.v1 < v1)
    return true;
  else if (aPoint.v1 > v1)
    return false;
  else if (aPoint.v1 == v1) {
    // start comparing v2
    if (aPoint.v2 < v2)
      return true;
    else if (aPoint.v2 > v2)
      return false;
    else if (aPoint.v2 == v2) {
      // compare v3
      if (aPoint.v3 < v3)
	return true;
      else if (aPoint.v3 >= v3) // notice greater or equal here
	return false;
    }
  }
}
