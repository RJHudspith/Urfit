/**
   @file Mainfile.c
   @brief Mainly a mainfile
 */
#include "gens.h"

#include "an_wrapper.h"
#include "io_wrapper.h"

#include "init.h"
#include "read_inputs.h"

#include "autocorr.h"
#include "bin.h"
#include "stats.h"

#include "write_flat.h"

//#define VIEW_AUTOCORR

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

  // bin the data we have read in
  if( bin_data( &Input ) == FAILURE ) {
    goto free_failure ;
  }

  // have a look at the time series if we have FFTW3
#ifdef VIEW_AUTOCORR
  #ifdef HAVE_FFTW3_H
  if( ACmeasure( Input ) == FAILURE ) {
    goto free_failure ;
  }
  #endif
#endif

  /*
  if( write_flat_file( Input , "Pseudo_LL_EM" ) == FAILURE ) {
    goto free_failure ;
  }
  */

  // resample the data we have read in
  if( resample_data( &Input ) == FAILURE ) {
    goto free_failure ;
  }

  /*
  if( write_flat_file( Input , "VVTitTit_EM_lc.flat" ) == FAILURE ) {
    goto free_failure ;
  }
  */

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
