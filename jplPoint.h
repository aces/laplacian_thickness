/* class jplPoint:
   Just a test class for now. Ultimately, if I stick with this scheme, I'll
   turn this into a template class

   Right now it just has three public int members
*/


#ifndef _JPL_POINT_H_
#define _JPL_POINT_H_

#include <iostream>

class jplPoint {
private:
  // nothing here yet
public:
  short v1;
  short v2;
  short v3;
  jplPoint(); //constructor
  jplPoint(short point1, short point2, short point3); // constructor
  void setPoints(short point1, short point2, short point3);
  bool operator<(const jplPoint& aPoint) const;
  friend ostream& operator<< (ostream& out, const jplPoint& aPoint);
};

#endif // _JPL_POINT_H_

  
