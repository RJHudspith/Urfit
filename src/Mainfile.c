/**
   @file Mainfile.c
   @brief Mainly a mainfile
 */
#include "gens.h"

#include "fake.h"
#include "fit_and_plot.h"
#include "init.h"
#include "pmap.h"
#include "read_corr.h"
#include "read_inputs.h"

#include "bootstrap.h"
#include "jacknife.h"
#include "rng.h"
#include "stats.h"

int
main( const int argc , const char *argv[] )
{
  if( argc != 3 ) {
    printf( "USAGE :: ./URFIT -i infile" ) ;
    return -1 ;
  }
  size_t i ;

  struct input_params Input ;
  if( read_inputs( &Input , argv[2] ) == FAILURE ) {
    goto free_failure ;
  }

  switch( Input.FileType ) {
  case Corr :
    if( read_corr( &Input , NOFOLD , 5 , 5 ) == FAILURE ) {
      goto free_failure ;
    } break ;
    // set Lt
    if( init_LT( &Input.Data , Input.Traj ) == FAILURE ) {
      goto free_failure ;
    }
  case Fake :
    if( generate_fake_data( &Input.Data , Input.Fit ,
			    Input.Traj , 0.0 , 0.005 ) == FAILURE ) {
      goto free_failure ;
    } break ;
  }
  
  // bootstrap it
  switch( Input.Data.Restype ) {
  case BootStrap :
    init_rng( 123456 ) ;
    bootstrap_full( &Input ) ;
    free_rng() ;
    break ;
  case JackKnife :
    jackknife_full( &Input ) ;
    break ;
  case Raw :
    for( i = 0 ; i < Input.Data.Ntot ; i++ ) {
      compute_err( &Input.Data.x[i] ) ;
      compute_err( &Input.Data.y[i] ) ;
    }
    break ;
  }
  
  // need to set this after data has been read ...
  if( fit_and_plot( Input ) == FAILURE ) {
    goto free_failure ;
  }

 free_failure :

  // free the structs
  free_Data( &Input.Data , Input.Fit ) ;
  free_Fit( &Input.Fit , Input.Data ) ;
  
  free_inputs( &Input ) ;
  
  return SUCCESS ;
}
