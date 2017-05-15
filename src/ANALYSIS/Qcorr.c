/**
   @file Qcorr.c
   @brief topological correlator measurement
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"

int
fit_Qcorr( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    const double OneV = 1.0 / ( Input -> Traj[i].Dimensions[0] * 
				Input -> Traj[i].Dimensions[1] *
				Input -> Traj[i].Dimensions[2] *
				Input -> Traj[i].Dimensions[3] ) ;

    const double t02 = 4.0 * ( Input -> Traj[i].Dimensions[0] *
			       Input -> Traj[i].Dimensions[0] ) / 100. ;

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      /*
      // flip the sign
      mult_constant( &Input -> Data.y[i] , -1 ) ;

      // compute r ... x is r^2
      root( &Input -> Data.x[i] ) ;
    
      // multiply by r
      mult( &Input -> Data.y[i] , Input -> Data.x[i] ) ;
      */
      mult_constant( &Input -> Data.y[j] , t02 * OneV ) ;
    }
    shift = j ;
  }
  
  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
	
  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}

int
fit_Qsusc( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // compute the partial sum up to whatever cutoff
    for( j = shift+1 ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      
      add( &Input -> Data.y[j] , Input -> Data.y[j-1] ) ;

    }
    
    shift = j ;
  }

  // normalise
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim; i++ ) {
    
    const double One_V = 1.0 / ( Input -> Traj[i].Dimensions[0] *
				 Input -> Traj[i].Dimensions[1] *
				 Input -> Traj[i].Dimensions[2] *
				 Input -> Traj[i].Dimensions[3] ) ;

    const double t02 = 4. * ( Input -> Traj[i].Dimensions[0] *
			      Input -> Traj[i].Dimensions[0] ) / 100. ;
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , One_V * t02 ) ;
    }
    shift = j ;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
	
  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}
