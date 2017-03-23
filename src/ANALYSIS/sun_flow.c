/**
   @file sun_flow.c
   @brief gradient flow analysis
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "resampled_ops.h"

int
sun_wflow_analysis( struct input_params *Input )
{
  // shift the data x = x - x[0]
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t idx = 0 , LT[9] = { 10,11,12,13,14,15,15,16,0} ;
    struct resampled tmp = init_dist( &(Input -> Data.x[shift]) ,
				      Input -> Data.x[shift].NSAMPLES , 
				      Input -> Data.x[shift].restype ) ;
    equate( &tmp , Input -> Data.x[shift] ) ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      subtract( &(Input -> Data.x[j]) , tmp ) ;
    }

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      divide_constant( &(Input -> Data.y[j]) , (double)LT[idx] ) ;
      idx++ ;
    }
    
    shift = j ;
    free( tmp.resampled ) ;
  }

  // perform the fit
  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  
  // free the fit
  if( fit != NULL ) {
    for( i = 0 ; i < Input -> Fit.Nparam ; i++ ) {
      if( fit[i].resampled != NULL ) {
	free( fit[i].resampled ) ;
      }
    }
    free( fit ) ;
  }
  
  return SUCCESS ;
}
