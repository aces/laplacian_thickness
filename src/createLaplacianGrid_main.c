#include "createLaplacianGrid.h"
#include  <ParseArgv.h>

/* command-line argument default values */
int     include_white_boundary     = 0;
int     include_grey_boundary       = 0;

/* argument parsing table */
ArgvInfo argTable[] = {
  { "-include_white_boundary", ARGV_CONSTANT, (char *)1, 
    (char *) &include_white_boundary,
    "Force the inclusion of the white-matter surface boundary. [Default: false]" },
  { "-include_grey_boundary", ARGV_CONSTANT, (char *)1,
    (char *) &include_grey_boundary,
    "Force the inclusion of the grey-matter surface boundary. [Default: false]"},
  
  { NULL, ARGV_END, NULL, NULL, NULL }
};

int main( int argc, char *argv[] ) {

  STRING         input_volume_filename, output_file_name;
  STRING         grey_surface_filename, white_surface_filename;
  STRING         history;

  Volume         out_volume;
  int            sizes[MAX_DIMENSIONS];
  int            x,y,z,i;

  /* create the history string from the input arguments */
  /* should check for overflow */
  history = alloc_string(1024); 
  for (i=0; i<argc; i++) {
    strcat(history, argv[i]);
    strcat(history, " ");
  }
  
  if ( ParseArgv( &argc, argv, argTable, 0 ) || ( argc != 5 )) {
    print_error("Usage: %s in_volume grey_surface.obj white_surface.obj output.mnc \n", argv[0] );
    return( 1 );
  }
  
  input_volume_filename  = argv[1];
  grey_surface_filename  = argv[2];
  white_surface_filename = argv[3];
  output_file_name       = argv[4];
  

  if( create_mantle( input_volume_filename, grey_surface_filename,
                     white_surface_filename, include_white_boundary,
                     include_grey_boundary, out_volume ) != OK ) {
    print_error( "Error: Cannot create mantle %s for Laplacian field\n", output_file_name );
  } else {
    output_volume(output_file_name, NC_BYTE, FALSE, 0.0, 0.0,
		  out_volume, history, NULL );
    printf( "Created mantle %s for Laplacian field\n", output_file_name );
  }

}
