/**
   @file init.c
   @brief initialise and free our fit and data structs
 */
#include "gens.h"

#include "pmap.h"

void
free_Data( struct data_info *Data ,
	   struct fit_info Fit )
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
    if( Fit.Corrfit == CORRELATED ) {
      if( Data -> Cov.W != NULL ) {
	if( Data -> Cov.W[i] != NULL ) {
	  free( Data -> Cov.W[i] ) ;
	}
      }
    }
  }
  if( Fit.Corrfit == UNCORRELATED ) {
    if( Data -> Cov.W != NULL ) {
      if( Data -> Cov.W[0] != NULL ) {
	free( Data -> Cov.W[0] ) ;
      }
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
  if( Fit -> Sims != NULL ) {
    free( Fit -> Sims ) ;
  }
  
  // free the priors
  if( Fit -> Prior != NULL ) {
    free( Fit -> Prior ) ;
  }
  return ;
}
