#include "corticalMantle.h"

#include <iostream.h>

void main() {
  corticalMantle *testMantle = new corticalMantle("grey.obj", "white.obj", "out.mnc", 
                                  "cls.mnc");
  testMantle->scanObjectsToVolume();
  testMantle->outputVolume("testout.mnc");
}

