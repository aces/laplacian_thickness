#include <iostream.h>
#include "corticalMantle.h"

extern "C" {
#include "ParseArgv.h"
}

// set up argument parsing defaults
int outsideValue = 10000;
int insideValue = 0;
int mantleValue = 5000;

// argument parsing table
ArgvInfo argTable[] = {
  { "-outside", ARGV_INT, (char *)0, (char*)&outsideValue,
    "Value to use outside of cortex" },
  { "-inside", ARGV_INT, (char*)0, (char*)&insideValue,
    "Value to use inside of cortex" },
  { "-mantle", ARGV_INT, (char*)0, (char*)&mantleValue,
    "Value to use for the mantle" },

  { NULL, ARGV_END, NULL, NULL, NULL }
};

//! Extracts the cortical mantle and initialises for thickness solving

int main(int argc, char *argv[]) {
  corticalMantle *mantle;

  if ( ParseArgv(&argc, argv, argTable, 0) || ( argc != 5 )) {
    cerr << "Usage: create_mantle [options] outer_surface.obj inner_surface.obj output.mnc like_this_file.mnc" << endl;
    exit(1);
  }
  
  char* outerSurfaceFile = argv[1];
  char* innerSurfaceFile = argv[2];
  char* outputFile = argv[3];
  char* likeFile = argv[4];

  mantle = new corticalMantle(outerSurfaceFile, innerSurfaceFile,
                              outputFile, likeFile);

  mantle->setVerbosity(1); // should be made an argument?
  mantle->scanObjectsToVolume();

  //  mantle->initialiseLaplacianGrid(outsideValue, insideValue, mantleValue,
  //FALSE, TRUE);

  short whiteInts, greyInts;

  mantle->countIntersections(76,120,71,greyInts, whiteInts);
  cout << "White: " << whiteInts << " Grey: " << greyInts << endl;

  mantle->output(outputFile);
}

