#include  <internal_volume_io.h>
#include  <bicpl.h>

int  main(
    int   argc,
    char  *argv[] )
{
    STRING               input_volume_filename;
    STRING               grey_object_filename, white_object_filename;
    STRING               output_filename;
    File_formats         grey_format, white_format;
    Volume               volume;
    int                  point, grey_n_points, grey_n_objects;
    int                  white_n_points, white_n_objects;
    Point                *grey_points, *white_points;
    object_struct        **grey_objects, **white_objects;
    Real                 value;
    FILE                 *file;

    initialize_argument_processing( argc, argv );

    if( !get_string_argument( NULL, &input_volume_filename ) ||
        !get_string_argument( NULL, &grey_object_filename ) ||
        !get_string_argument( NULL, &white_object_filename ) ||
        !get_string_argument( NULL, &output_filename ) )
    {
        print_error(
           "Usage: %s  volume.mnc grey_object.obj white_object.obj output_file.txt\n", argv[0]);
        return( 1 );
    }

    if( input_volume( input_volume_filename, 3, File_order_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &volume, NULL ) != OK )
        return( 1 );

    if( input_graphics_file( grey_object_filename,
                             &grey_format, &grey_n_objects, 
                             &grey_objects ) != OK )
        return( 1 );

    grey_n_points = get_object_points( grey_objects[0], &grey_points );


    if( input_graphics_file( white_object_filename,
                             &white_format, &white_n_objects, 
                             &white_objects ) != OK )
        return( 1 );

    white_n_points = get_object_points( white_objects[0], &white_points );


    if( open_file( output_filename, WRITE_FILE, ASCII_FORMAT, &file ) != OK )
        return( 1 );

    for_less( point, 0, white_n_points )
    {


        evaluate_volume_in_world( volume,
                                  (RPoint_x(white_points[point]) +
                                   RPoint_x(grey_points[point])) / 2,
                                  (RPoint_y(white_points[point]) +
                                   RPoint_y(grey_points[point])) / 2,
                                  (RPoint_z(white_points[point]) +
                                   RPoint_z(grey_points[point])) / 2,
                                  0, FALSE, 0.0,
                                  &value,
                                  NULL, NULL, NULL,
                                  NULL, NULL, NULL, NULL, NULL, NULL );

        if( output_real( file, value ) != OK ||
            output_newline( file ) != OK )
            return( 1 );
    }

    (void) close_file( file );

    delete_object_list( white_n_objects, white_objects );
    delete_object_list( grey_n_objects, grey_objects );

    return( 0 );
}
