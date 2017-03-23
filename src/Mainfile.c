/**
   @file Mainfile.c
   @brief Mainly a mainfile
 */
#include "gens.h"

#include "an_wrapper.h"
#include "io_wrapper.h"

#include "fit_and_plot.h"
#include "init.h"

#include "read_inputs.h"


#include "bootstrap.h"
#include "jacknife.h"
#include "rng.h"
#include "stats.h"

int
main( const int argc , const char *argv[] )
{
  // tell us if we did something stupid
  if( argc != 3 ) {
    fprintf( stderr , "USAGE :: ./URFIT -i infile" ) ;
    return -1 ;
  }
  size_t i ;

  // initially read the inputs
  struct input_params Input ;
  if( read_inputs( &Input , argv[2] ) == FAILURE ) {
    goto free_failure ;
  }

  // choose the correct IO
  if( io_wrap( &Input ) == FAILURE ) {
    goto free_failure ;
  }
  
  // bootstrap it
  switch( Input.Data.Restype ) {
  case BootStrap :
    printf( "Bootstrapping \n" ) ;
    init_rng( 123456 ) ;
    bootstrap_full( &Input ) ;
    free_rng() ;
    break ;
  case JackKnife :
    printf( "JackKnifing \n" ) ;
    jackknife_full( &Input ) ;
    break ;
  case Raw :
    printf( "Raw data\n" ) ;
    for( i = 0 ; i < Input.Data.Ntot ; i++ ) {
      compute_err( &Input.Data.x[i] ) ;
      compute_err( &Input.Data.y[i] ) ;
    }
    break ;
  }

  printf( "Bootstrapping finished\n" ) ;
  
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
