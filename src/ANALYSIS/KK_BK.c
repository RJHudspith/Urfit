/**
   @file KKops.c
   @brief code computes the ratio O_i / O_1
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "resampled_ops.h"
#include "init.h"
#include "effmass.h"

#define SYMMETRIZE1

#define SYMMETRIZE2

int
fit_bags( struct input_params *Input )
{  
  size_t i , j ;
  const size_t LT = Input -> Traj[0].Dimensions[3] ;
  const size_t SRC = Input -> Traj[0].Dimensions[0] ;

  const double V = Input -> Traj[0].Dimensions[1] *
    Input -> Traj[0].Dimensions[1] *
    Input -> Traj[0].Dimensions[1] ;

  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    // compute the denominator
    mult( &Input -> Data.y[j+LT] , Input -> Data.y[(LT-j+SRC)%LT+2*LT] ) ;
    divide( &Input -> Data.y[j] , Input -> Data.y[j+LT] ) ;
    mult_constant( &Input -> Data.y[j] , -2*V ) ;
    //equate( &Input -> Data.y[j] , Input -> Data.y[j+LT] ) ;
  }

  #ifdef SYMMETRIZE1
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    if( j < LT/2 ) {
      add( &Input -> Data.y[j] , Input -> Data.y[j + LT/2 ] ) ;
      mult_constant( &Input -> Data.y[j] , 0.5 ) ;
    } else {
      mult_constant( &Input -> Data.y[j] , 0.0 ) ;
    }
  }
    #ifdef SYMMETRIZE2
    for( j = 1 ; j < Input -> Data.Ndata[0] ; j++ ) {
      if( j <= SRC/2 ) {
	add( &Input -> Data.y[j] , Input -> Data.y[(SRC-j)%SRC ] ) ;
	mult_constant( &Input -> Data.y[j] , 0.5 ) ;
      } else {
	mult_constant( &Input -> Data.y[j] , 0.0 ) ;
      }
    }
    #endif
  #endif

  Input -> Data.Nsim = Input -> Data.Nsim/3 ;
  Input -> Data.Ntot = Input -> Data.Ntot/3 ;
  
  struct resampled *effmass = effective_mass( Input , LOGBWD_EFFMASS ) ;
  
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;
  
  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;  

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  Input -> Data.Nsim = Input -> Data.Nsim*3 ;
  Input -> Data.Ntot = Input -> Data.Ntot*3 ;
  
  return SUCCESS ;
}

