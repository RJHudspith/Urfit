/**
   @file init.c
   @brief initialise and free our fit and data structs
 */
#include "gens.h"

void
free_Data( struct data_info *Data )
{
  // free all the data
  size_t i ;
  for( i = 0 ; i < Data -> Ntot ; i++ ) {
    if( Data -> x[i].resampled != NULL ) {
      free( Data -> x[i].resampled ) ;
    }
    if( Data -> y[i].resampled != NULL ) {
      free( Data -> y[i].resampled ) ;
    }
    if( Data -> Cov.W[i] != NULL ) {
      free( Data -> Cov.W[i] ) ;
    }
  }
  if( Data -> x != NULL ) {
    free( Data -> x ) ;
  }
  if( Data -> y != NULL ) {
    free( Data -> y ) ;
  }
  if( Data -> Cov.W != NULL ) {
    free( Data -> Cov.W ) ;
  }
  if( Data -> Ndata != NULL ) {
    free( Data -> Ndata ) ;
  }
  return ;
}

void
free_Fit( struct fit_info *Fit ,
	  const struct data_info Data )
{
  // free the simultaneous parameter map
  free( Fit -> Sims ) ;

  // free the full parameter map
  free_pmap( Fit -> map , Data.Ntot ) ;
}