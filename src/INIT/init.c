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
  }
  if( Data -> Cov.W != NULL ) {
    if( Fit.Corrfit == UNCORRELATED ) {
      if( Data -> Cov.W[0] != NULL ) {
	free( Data -> Cov.W[0] ) ;
      }
    } else if( Fit.Corrfit == UNCORRELATED ) {
      for( i = 0 ; i < Data -> Ntot ; i++ ) {
	if( Data -> Cov.W[i] != NULL ) {
	  free( Data -> Cov.W[i] ) ;
	}
      }
    }
    free( Data -> Cov.W ) ;
  }
  if( Data -> x != NULL ) {
    free( Data -> x ) ;
  }
  if( Data -> y != NULL ) {
    free( Data -> y ) ;
  }
  if( Data -> Ndata != NULL ) {
    free( Data -> Ndata ) ;
  }
  if( Data -> LT != NULL ) {
    free( Data -> LT ) ;
  }
  return ;
}

// free the fit information struct
void
free_Fit( struct fit_info *Fit )
{
  // free the simultaneous parameter map
  if( Fit -> Sims != NULL ) {
    free( Fit -> Sims ) ;
  }
  
  // free the priors
  if( Fit -> Prior != NULL ) {
    free( Fit -> Prior ) ;
  }
  // free the guesses
  if( Fit -> Guess != NULL ) {
    free( Fit -> Guess ) ;
  }
  return ;
}

// free the fit parameter array
void
free_fitparams( struct resampled *Fit ,
		const size_t Nlogic )
{
  size_t i ;
  if( Fit != NULL ) {
    for( i = 0 ; i < Nlogic ; i++ ) {
      if( Fit[i].resampled != NULL ) {
	free( Fit[i].resampled ) ;
      }
    }
    free( Fit ) ;
  }
  return ;
}

// initialise the "LT" array
int
init_LT( struct data_info *Data ,
	const struct traj *Traj )
{
  size_t shift = 0 , j , i ;
  Data -> LT = malloc( Data -> Ntot * sizeof( size_t ) ) ;
  for( i = 0 ; i < Data -> Nsim ; i++ ) {
    for( j = shift ; j < shift + Data -> Ndata[i] ; j++ ) {
      Data -> LT[j] = Traj[i].Dimensions[ Traj[i].Nd - 1 ] ;
    }
    shift = j ;
  }
  return SUCCESS ;
}
