/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/
/*
 * Code adapted from David MacDonald's surface_mask2
 */

#include "createLaplacianGrid.h"


#define  BINTREE_FACTOR   0.1

#define  OFFSET   1000.0

private  BOOLEAN  get_extrapolated_value(
    Volume         volume,
    int            x,
    int            y,
    int            z,
    object_struct  *surface,
    Real           threshold,
    Real           *value_to_set );


int binary_object_mask( STRING input_surface_filename,
                        Real set_value,
                        Real inside_value,
                        Volume volume ) 
{

    Real                 separations[MAX_DIMENSIONS];
    Real                 value;
    Real                 xw, yw, zw, value_to_set;
    Real                 sign_before, sign_after;
    Point                origin, dest;
    Real                 min_real_value, max_real_value;
    Real                 *directions;
    Vector               direction;
    STRING               history;
    File_formats         format;
    int                  x, y, z, n_objects, k;
    int                  sizes[MAX_DIMENSIONS];
    int                  obj_index;
    int                  i, j, best, n_intersections, ind;
    Real                 dist, *distances, tmp, voxel[MAX_DIMENSIONS];
    object_struct        **objects;
    BOOLEAN              set_value_specified, inside, erase_flag;
    BOOLEAN              threshold_specified;
    progress_struct      progress;


    if( input_graphics_file( input_surface_filename,
                             &format, &n_objects, &objects ) != OK )
        return( 1 );

    if( n_objects < 1 )
    {
        print( "No objects in %s.\n", input_surface_filename);
        return( 1 );
    }

    get_volume_separations( volume, separations );
    get_volume_sizes( volume, sizes );

    min_real_value = get_volume_real_min( volume );
    max_real_value = get_volume_real_max( volume );

    if( BINTREE_FACTOR > 0.0 )
    {
        create_polygons_bintree( get_polygons_ptr(objects[0]),
                     ROUND( (Real) get_polygons_ptr(objects[0])->n_items *
                                                   BINTREE_FACTOR ) );
    }

    initialize_progress_report( &progress, FALSE, sizes[X] * sizes[Y],
                                "Masking" );

    for_less( x, 0, sizes[X] )
    {
        for_less( y, 0, sizes[Y] )
        {
            voxel[X] = (Real) x;
            voxel[Y] = (Real) y;
            voxel[Z] = -OFFSET;

            convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );
            fill_Point( origin, xw, yw, zw );

            voxel[Z] = (Real) sizes[Z] + OFFSET;
            convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );
            fill_Point( dest, xw, yw, zw );
            SUB_POINTS( direction, dest, origin );
            NORMALIZE_VECTOR( direction, direction );

            initialize_intersect_directions();

            n_intersections = intersect_ray_with_object( &origin, &direction,
                                       objects[0], &obj_index, &dist,
                                       &distances );

            directions = get_intersect_directions();

            for_less( i, 0, n_intersections )
                distances[i] = distances[i] / FABS( separations[Z] ) - OFFSET;

            for_less( i, 0, n_intersections-1 )
            {
                best = i;
                for_less( j, i+1, n_intersections )
                {
                    if( distances[j] < distances[best] )
                        best = j;
                }

                tmp = distances[i];
                distances[i] = distances[best];
                distances[best] = tmp;
                tmp = directions[i];
                directions[i] = directions[best];
                directions[best] = tmp;
            }

            j = 0;
            i = 0;
            while( i < n_intersections )
            {
                if( j == 0 || distances[i] - distances[j-1] > 1e-10 ||
                    directions[i] != directions[j-1] )
                {
                    distances[j] = distances[i];
                    directions[j] = directions[i];
                    ++j;
                }
                ++i;
            }
            n_intersections = j;

            while( TRUE )
            {
                for_less( i, 0, n_intersections-1 )
                {
                    if( distances[i+1] - distances[i] < 1e-10 )
                        break;
                }

                if( i >= n_intersections-1 )
                    break;

                j = i+1;
                while( j < n_intersections-1 &&
                       distances[j+1] - distances[i] < 1e-10 )
                    ++j;

                if( i == 0 )
                    sign_before = 1.0;
                else
                    sign_before = directions[i-1];

                if( j == n_intersections-1 )
                    sign_after = -1.0;
                else
                    sign_after = directions[j+1];

                if( sign_before == sign_after )
                {
                    for_less( k, j, n_intersections )
                    {
                        distances[i+k-j] = distances[k];
                        directions[i+k-j] = directions[k];
                    }
                    directions[j] = -sign_before;
                    n_intersections -= j - i;
                }
                else
                {
                    for_less( k, j+1, n_intersections )
                    {
                        distances[i+k-j-1] = distances[k];
                        directions[i+k-j-1] = directions[k];
                    }
                    n_intersections -= j - i + 1;
                }
            }

            if( n_intersections % 2 == 0 )
            {
                ind = 0;
                inside = FALSE;
                for_less( z, 0, sizes[Z] )
                {
                    while( ind <= n_intersections-1 &&
                           (Real) z >= distances[ind] )
                    {
                        inside = !inside;
                        ++ind;
                    }

                    erase_flag = FALSE;
                    if( !inside )
                    {
                      erase_flag = TRUE;
                    }
                    if( erase_flag )
                    {
                      value_to_set = set_value;
                      set_volume_real_value( volume, x, y, z, 0, 0,
                                             value_to_set );
                    }
                    else {
                      set_volume_real_value( volume, x, y, z, 0, 0, 
                                             inside_value );
                    }
                }
            }
            else {
                print( "N intersections: %d\n", n_intersections );
            }

            if( n_intersections > 0 )
            {
                FREE( distances );
                FREE( directions );
            }

            update_progress_report( &progress, x * sizes[Y] + y + 1 );
        }
    }

    terminate_progress_report( &progress );
    
    delete_object_list( n_objects, objects );
    return( 0 );
}

private  BOOLEAN  get_one_extrapolated_value(
    Volume         volume,
    int            x1,
    int            y1,
    int            z1,
    int            x2,
    int            y2,
    int            z2,
    object_struct  *surface,
    Real           threshold,
    Real           *value_to_set )
{
    int      sizes[N_DIMENSIONS];
    int      obj_index, n_int, n_intersections, i;
    Real     value2, dist, *distances, xw, yw, zw;
    Real     voxel[N_DIMENSIONS];
    Point    p1, p2;
    Vector   direction;

    get_volume_sizes( volume, sizes );

    if( x2 < 0 || x1 >= sizes[X] ||
        y2 < 0 || y1 >= sizes[Y] ||
        z2 < 0 || z1 >= sizes[Z] )
        return( FALSE );

    voxel[X] = (Real) x1;
    voxel[Y] = (Real) y1;
    voxel[Z] = (Real) z1;
    convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );
    fill_Point( p1, xw, yw, zw );

    voxel[X] = (Real) x2;
    voxel[Y] = (Real) y2;
    voxel[Z] = (Real) z2;
    convert_voxel_to_world( volume, voxel, &xw, &yw, &zw );
    fill_Point( p2, xw, yw, zw );
    SUB_POINTS( direction, p2, p1 );

    n_intersections = intersect_ray_with_object( &p1, &direction,
                                                 surface, &obj_index, &dist,
                                                 &distances );

    n_int = 0;
    for_less( i, 0, n_intersections )
    {
        if( distances[i] > 0.0 && distances[i] < 1.0 )
        {
            ++n_int;
            dist = distances[i];
        }
    }

    if( n_intersections > 0 )
        FREE( distances );

    if( n_int != 1 )
        return( FALSE );

    value2 = get_volume_real_value( volume, x2, y2, z2, 0, 0 );

    if( value2 <= threshold )
        return( FALSE );

    *value_to_set = (threshold - value2 * dist) / (1.0 - dist);

    return( TRUE );
}

private  BOOLEAN  get_extrapolated_value(
    Volume         volume,
    int            x,
    int            y,
    int            z,
    object_struct  *surface,
    Real           threshold,
    Real           *value_to_set )
{
    Real     value;
    int      c, dir, n_bounds, voxel[N_DIMENSIONS];

    *value_to_set = 0.0;
    n_bounds = 0;
    for_less( c, 0, N_DIMENSIONS )
    {
        for( dir = -1;  dir <= 1;  dir += 2 )
        {
            voxel[X] = x;
            voxel[Y] = y;
            voxel[Z] = z;
            voxel[c] += dir;
            if( get_one_extrapolated_value( volume, x, y, z,
                                            voxel[X], voxel[Y], voxel[Z],
                                            surface, threshold, &value ) )
            {
                *value_to_set += value;
                ++n_bounds;
            }
        }
    }

    if( n_bounds > 0 )
        *value_to_set /= (Real) n_bounds;

    return( n_bounds > 0 );
}

int create_mantle ( char *input_volume_filename,
		    char *grey_surface_filename,
		    char *white_surface_filename,
		    int include_white_boundary,
		    int include_grey_boundary,
                    Volume * out_volume )
{
  /*
    STRING         input_volume_filename, output_file_name;
    STRING         grey_surface_filename, white_surface_filename;
    STRING         history;
  */
  Volume         grey_volume , white_volume;
  int            sizes[MAX_DIMENSIONS];
  int            x,y,z;
  Real           white_value, grey_value;
  object_struct  **objects;  
  File_formats   format;
  int            n_objects, i;
  Real outside_value = 0;
  Real inside_value = 1;

  int ndims;
  STRING * dimension_names;
  get_file_dimension_names( input_volume_filename, &ndims, &dimension_names );
  if( ndims != 3 ) {
    fprintf( stderr, "Error: Dimension of input volume %s must be 3.\n", 
             input_volume_filename );
    return( 1 );
  }

  if( input_volume( input_volume_filename, 3, dimension_names,
                    NC_BYTE, FALSE, 0.0, 255.0, TRUE, 
                    &white_volume, NULL ) != OK ) {
    return( 1 );
  }
  set_volume_real_range(white_volume, 0, 10); 

  grey_volume = copy_volume_definition( white_volume, NC_BYTE,
                                        FALSE, 0.0, 255.0 );
  set_volume_real_range(grey_volume, 0, 10); 

  *out_volume = copy_volume_definition( white_volume, NC_BYTE,
                                        FALSE, 0.0, 255.0 );
  set_volume_real_range(*out_volume, 0, 10); 

  (void) binary_object_mask( grey_surface_filename,
                             outside_value,
                             inside_value,
                             grey_volume );

  (void) binary_object_mask( white_surface_filename,
                             outside_value,
                             inside_value,
                             white_volume );

  get_volume_sizes( *out_volume, sizes );

  /* now test for various possible conditions */
  for_less( x, 0, sizes[0] ) {
    for_less( y, 0, sizes[1] ) {
      for_less( z, 0, sizes[2] ) {
        grey_value = get_volume_real_value(grey_volume, x, y, z, 0, 0);
        white_value = get_volume_real_value(white_volume, x, y, z, 0, 0);
        /*        printf("Grey: %f White: %f\n", grey_value, white_value);*/
        if (grey_value > 0.5 && white_value > 0.5 ) {
          /* inside both surfaces */
          set_volume_real_value(*out_volume, x, y, z, 0, 0, 0);
        }
        else if (grey_value > 0.5 && white_value < 0.5) {
          /* in mantle */
          set_volume_real_value(*out_volume, x, y, z, 0, 0, 5);
        }
        else if (grey_value < 0.5 && white_value < 0.5) {
          /* outside both surfaces */
          set_volume_real_value(*out_volume, x, y, z, 0, 0, 10);
        }
        else if (grey_value < 0.5 && white_value > 0.5){
          /* shouldn't happen, as this means overlap - treat as mantle
             for now */
          set_volume_real_value(*out_volume, x, y, z, 0, 0, 2);
        }
      }
    }
  }

  delete_volume( grey_volume );
  delete_volume( white_volume );

  /* force inclusion of the boundaries if so desired */
  if ( include_white_boundary == 1 ) {
    input_graphics_file(white_surface_filename, &format, &n_objects,
                        &objects);
    scan_object_to_volume(objects[0], *out_volume, *out_volume, 5000, 1);
  }
  if ( include_grey_boundary == 1 ) {
    input_graphics_file(grey_surface_filename, &format, &n_objects,
                        &objects);
    scan_object_to_volume(objects[0], *out_volume, *out_volume, 5000, 1);
  }

  return OK;
}


