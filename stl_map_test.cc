/* Test code for me to get used to the STL's map template */


#include <iostream.h>
#include <map>
#include <list>
#include "jplPoint.h"

void main() {
  typedef map<jplPoint, int>::const_iterator point_iterator;

  list<jplPoint> testList(0);

  testList.push_back(jplPoint (1,2,3));
  testList.push_back(jplPoint (3,4,5));

  jplPoint testG = testList.front();
  cout << "Vector Point: " << testG << endl;
  testList.pop_front();
  jplPoint testF = testList.front();
  cout << "Vector Point: " << testF << endl;

  jplPoint test(3, 8, 5);
  jplPoint test2(2, 1, 3);

  cout << "Middle Point: " << test.v2 << endl;

  map<jplPoint, int> testMap;

  testMap[test] = 5;
  testMap[test2] = 2;
  testMap[jplPoint (2,2,2)] = 1;

  cout << "First test case: " << testMap[jplPoint (3,8,5)] << endl;

  cout << "Size: " << testMap.size() << endl;

  for (point_iterator m = testMap.begin(); m != testMap.end(); ++m) {
    jplPoint test5 = m->first;
    cout << "First v: " << test5.v1;
    cout << "Value: " << m->first << " Second: " << m->second <<  endl;
  }

}

  
