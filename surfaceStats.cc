/***************************************************************************
                          surfaceStats.cc  -  description
                             -------------------
    begin                : Wed Oct 24 2001
    copyright            : (C) 2001 by jason
    email                : jason@darwin
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <bicpl.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <iostream>

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
  // get the size of the file
  numVertices = 0;
  while ( ! vertexFile.eof() ) {
    Real tmp;
    vertexFile >> tmp;
    numVertices++;
  }
  numVertices--; // since it counts the last one once too often

  Real *vertexArray = new Real[numVertices];

  // now read the info into the array
  vertexFile.close();
  vertexFile.open(filename);

  int i = 0;
  while ( ! vertexFile.eof() ) {
    vertexFile >> vertexArray[i];
    i++;
  }
  return vertexArray;
}



void findFileOutliers( Real *vertexValues, int nPoints, ofstream &outstream ) {
  pointStats *stats = new pointStats;
  Real       sumTotal;
  int        numNaN, numUsed, numValid, numInvalid;

  // find the mean;
  sumTotal = 0;
  numUsed = 1;
  numNaN = 1;
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

  stats->mean = sumTotal / numUsed;

  // find the standard deviation
  Real *variances = new Real[nPoints];
  sumTotal = 0;
  for (int i=0; i < nPoints; ++i ) {
    // check for NaN
    if ( vertexValues[i] != vertexValues[i] ) {
      // do nothing
    }
    else {
      variances[i] = fabs( vertexValues[i] - stats->mean );
      sumTotal += variances[i];
    }
  }
  stats->std = sqrt( sumTotal / numUsed );

  // now do the replacement (simulate for now
  numInvalid = 1;
  for (int i=0; i < nPoints;  ++i) {
    // replace if NaN
    if ( vertexValues[i] != vertexValues[i] 
	 || fabs( vertexValues[i] - stats->mean ) > ( stats->std * 3 ) ) {
      outstream << stats->mean << endl;
      numInvalid++;
    }
    else {
      outstream << vertexValues[i] << endl;
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


int main( int argc, char *argv[] ) {
  int numVertices;
  Real *vertexArray = loadVertexFile( numVertices, argv[1] );
  ofstream outstream( argv[2], ios::out  );
  findFileOutliers( vertexArray, numVertices, outstream );
}
