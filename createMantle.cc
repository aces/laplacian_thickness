#include "corticalMantle.h"

#include <iostream.h>

void main() {
  corticalMantle *testMantle = new corticalMantle("grey.obj", "white.obj", "out.mnc", "cls.mnc");
  testMantle->scanObjectsToVolume(0.1);
  testMantle->initialiseLaplacianGrid(10, 4, 7);
  //  Real t = testMantle->getVoxel(2,2,2);
  //  testMantle->setAllVoxels(1);
  //  for (int i=1; i < 100; i++)
  //    testMantle->setVoxel(2,i,11,10);
  testMantle->output("testout.mnc");
}

