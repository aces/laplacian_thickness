#include "corticalMantle.h"

#include <iostream.h>

void main() {
  corticalMantle *testMantle = new corticalMantle("grey.obj", "white.obj", "out.mnc", 
                                  "cls.mnc");
  testMantle->scanObjectsToVolume();
  testMantle->initialiseLaplacianGrid(10, 4, 7);
  testMantle->outputVolume("testout.mnc");
}

