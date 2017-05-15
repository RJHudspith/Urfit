/**
   @file Mainfile.c
   @brief Mainly a mainfile
 */
#include "gens.h"

#include "an_wrapper.h"
#include "io_wrapper.h"

#include "init.h"
#include "read_inputs.h"

#include "histogram.h"

#include "stats.h"

int
main( const int argc , const char *argv[] )
{
  // tell us if we did something stupid
  if( argc != 3 ) {
    fprintf( stderr , "USAGE :: ./URFIT -i infile" ) ;
    return -1 ;
  }

  // initially read the inputs
  struct input_params Input ;
  if( read_inputs( &Input , argv[2] ) == FAILURE ) {
    goto free_failure ;
  }

  // choose the correct IO
  if( io_wrap( &Input ) == FAILURE ) {
    goto free_failure ;
  }

  // resample the data we have read in
  resample_data( &Input ) ;
  
  // need to set this after data has been read ...
  if( an_wrapper( &Input ) == FAILURE ) {
    goto free_failure ;
  }
  
 free_failure :

  // free the structs
  free_Data( &Input.Data , Input.Fit ) ;
  free_Fit( &Input.Fit , Input.Data ) ;
  
  free_inputs( &Input ) ;
  
  return SUCCESS ;
}
