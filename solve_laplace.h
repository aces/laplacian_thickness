#ifndef _SOLVE_LAPLACE_H_
#define _SOLVE_LAPLACE_H_

extern "C" {
#include "volume_io.h"
#include "bicpl.h"
}

#include <iostream.h>
#include <math.h>
#include "surface_volume_defs.h"

#define INNER_SURFACE 0
#define OUTER_SURFACE 10000
#define MANTLE_INIT 5000

void initialise_volume(Volume label_volume);
void solve_laplace(Volume laplacian_volume, int *sizes);
void gradient_volume(Volume laplace_volume, int *sizes);



#endif // _SOLVE_LAPLACE_H_
