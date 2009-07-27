/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/***************************************************************************
                          surfaceStats.cc  -  description
                             -------------------
    begin                : Wed Oct 24 2001
    email                : jason@darwin
 ***************************************************************************/

#include <bicpl.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <iostream>

extern "C" {
#include "ParseArgv.h"
}

using namespace std;

struct pointStats_struct {
  float mean;
  float std;
  int   numNeighbours;
};
typedef pointStats_struct pointStats;

struct neighbour_struct {
  int index;
  vector <int> neighbourIndices;
};
typedef neighbour_struct neighbour;

//neighbour *createNeighbourList( const polygon_struct &p ) {
//}

// load a vertex file, returning it in an 
Real *loadVertexFile( int &numVertices, char *filename ) {
  ifstream vertexFile(filename);
  string line;
  // get the size of the file
  numVertices = 0;
  while ( ! vertexFile.eof() ) {
    getline(vertexFile, line);
    numVertices++;
  }
  numVertices--; // since it counts the last one once too often

  Real *vertexArray = new Real[numVertices];

  // now read the info into the array
  vertexFile.close();
  vertexFile.open(filename);
  vertexFile.clear();
  
  for (int i=0; i< numVertices; i++) {
    getline(vertexFile, line);
    vertexArray[i] = atof(line.c_str());
  }
  return vertexArray;
}



void findFileOutliers( Real *vertexValues, int nPoints, 
		       bool replaceNaN=true,
		       bool replaceOutliers=true,
		       float outlierStdMultiplier=3.0) {
  pointStats *stats = new pointStats;
  Real       sumTotal;
  int        numNaN, numUsed, numValid, numInvalid;

  // find the mean;
  sumTotal = 0;
  numUsed = 0;
  numNaN = 0;
  for (int i=0; i < nPoints; ++i ) {
    // check for not a number
    if ( vertexValues[i] != vertexValues[i] ) {
      numNaN++;
    }
    else {
      sumTotal += vertexValues[i];
      numUsed++;
    }
  }

  if( numUsed > 0 ) {
    stats->mean = sumTotal / numUsed;
  } else {
    stats->mean = 0.0;
  }

  // find the standard deviation
  sumTotal = 0;
  for (int i=0; i < nPoints; ++i ) {
    // check for NaN
    if ( vertexValues[i] != vertexValues[i] ) {
      // do nothing
    }
    else {
      sumTotal += fabs( vertexValues[i] - stats->mean );
    }
  }
  if( numUsed > 1 ) {
    stats->std = sqrt( sumTotal / ( numUsed - 1 ) );
  } else {
    stats->std = 0.0;
  }

  // now do the replacement NOTE: instead of doing this globally it
  // ought to work based on neighbourhood information ...
  numInvalid = 0;
  for (int i=0; i < nPoints;  ++i) {
    // check if NaN
    if ( vertexValues[i] != vertexValues[i] ) {
      numInvalid++;
      if ( replaceNaN==true ) {
	vertexValues[i] = stats->mean;
      }
    }
    // check if value is an outlier, defined as the standard
    // deviation times supplied constant
    else if ( fabs( vertexValues[i] - stats->mean ) > 
	      ( stats->std * outlierStdMultiplier ) ) {
      numInvalid++;
      if ( replaceOutliers==true ) {
	vertexValues[i] = stats->mean;
      }
    }
  }

  cout << "STATS:" << endl 
       << "  # points = " << nPoints << endl 
       << "  # used   = " << numUsed << endl
       << "  # NaN    = " << numNaN << endl
       << "  # inval  = " << numInvalid << endl
       << "  mean     = " << stats->mean << endl
       << "  std      = " << stats->std << endl;
}

void writeVerticesToFile( char *filename, Real *vertexValues, int nPoints ) {
  ofstream outstream( filename, ios::out  );
  for ( int i=0; i < nPoints; ++i ) {
    outstream << vertexValues[i] << endl;
  }
}

// void findNeighbouringOutliers( const polygons_struct &p, Real *vertexValues ) {

//   for ( int f = 0; f < p.n_items; f++ ) {
//     pointStats *stats = new pointStats;
//     int numNeighbours = GET_OBJECT_SIZE( p, pointIndex );
//     Real *values = new Real[numNeighbours];
//     int *indices = new int[numNeighbours];
//     Real *variances = new Real[numNeighbours];
//     Real sumTotal = 0;

//     // generate an array containing all immediate neighbours
//     for ( int i = 0; i < numNeighbours; ++i ) {
//       indices[i] = p.indices[POINT_INDEX(p.end_indices, pointIndex, i)];
//       values[i] = vertexValues[indices[i]];
//       sumTotal += values[i];
//     }

//     // find the mean and standard deviation: still kinda silly since I'm for now
//     // only using immediate neighbours.
//     stats->mean = sumTotal / ( numNeighbours );
//     sumTotal = 0;
//     for ( i = 0; i < numNeighbours; ++i ) {
//       variances[i] = fabs( values[i] - stats->mean );
//       sumTotal += variances[i];
//     }
//     stats->std = sqrt( sumTotal );

//     for ( i = 0; i < numNeighbours; ++i ) {
//       if ( values[i] > ( stats->mean + (2 * stats->std) ) ||
//             values[i] < ( stats->mean - (2 * stats->std) ) ) {
//         cout << "Deviance at: " << indices[i] << "with a value of "
//               << values[i] << " a mean of " << stats->mean
//               << " and a std of " << stats->std << endl;
//       }
//     }
//   delete values;
//   delete indices;
//   delete variances;
//   }
// }
int replaceNaNs         = 0;
int replaceOutliers     = 0;
Real stdConstant        = 3.0;

ArgvInfo argTable[] = {
  { "-replace_nans", ARGV_CONSTANT, (char *) 1, (char *) &replaceNaNs,
    "Replace NaNs with mean." },
  { "-replace_outliers", ARGV_CONSTANT, (char *) 1,
      (char *) &replaceOutliers,
      "Replace outliers with mean." },
  { "-std_constant", ARGV_FLOAT, (char *)0, (char *) &stdConstant,
    "Define outlier as this constant times the standard deviation" },
  
  { NULL, ARGV_END, NULL, NULL, NULL }
};


int main( int argc, char *argv[] ) {

  if ( ParseArgv( &argc, argv, argTable, 0 ) || (! argc > 1 ) ) {
    cerr << "Usage: " << argv[0] << " [options] input.txt [output.txt]"
	 << endl;
    cerr << endl << "Copyright Alan C. Evans" << endl
                 << "Professor of Neurology" << endl
                 << "McGill University" << endl;
    return(1);
  }



  int numVertices;
  Real *vertexArray = loadVertexFile( numVertices, argv[1] );
  findFileOutliers( vertexArray, numVertices, replaceNaNs, replaceOutliers,
		    stdConstant);

  if ( replaceNaNs == TRUE || replaceOutliers == TRUE ) {
    if ( argc != 3 ) {
      cerr << "ERROR: output file missing" << endl;
      return(1);
    }
    writeVerticesToFile( argv[2], vertexArray, numVertices );
  }
}
