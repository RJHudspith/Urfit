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
#include "reweight.h"

#include "write_flat.h"

//#define VIEW_AUTOCORR
//#define VIEW_RAWDATA

int
main( const int argc , const char *argv[] )
{
  // tell us if we did something stupid
  if( argc != 3 ) {
    fprintf( stderr , "USAGE :: ./URFIT -i infile" ) ;
    return -1 ;
  }

  printf( "Here \n" ) ;

  // initially read the inputs
  struct input_params Input ;
  if( read_inputs( &Input , argv[2] ) == FAILURE ) {
    printf( "[INPUT] Input reading failed\n" ) ;
    goto free_failure ;
  }

  printf( "Input read \n" ) ;

  // choose the correct IO
  if( io_wrap( &Input ) == FAILURE ) {
    goto free_failure ;
  }

  // reweight the data
  if( reweight_data( &Input ) == FAILURE ) {
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

#ifdef VIEW_RAWDATA
  for( size_t i = 0 ; i < Input.Data.y[0].NSAMPLES ; i++ ) {
    for( size_t j = 0 ; j < Input.Data.Ndata[0] ; j++ ) {
      printf( "%e %e\n" ,
	      Input.Data.x[j].resampled[i] ,
	      Input.Data.y[j].resampled[i] ) ;
    }
    printf( "\n" ) ;
  }
#endif

  // resample the data we have read in
  if( resample_data( &Input ) == FAILURE ) {
    goto free_failure ;
  }

  // need to set this after data has been read ...
  if( an_wrapper( &Input ) == FAILURE ) {
    goto free_failure ;
  }
  
 free_failure :

  // free the structs
  free_Data( &Input.Data , Input.Fit ) ;
  free_Fit( &Input.Fit ) ;
  
  free_inputs( &Input ) ;
  
  return SUCCESS ;
}
