/**
   @file Qcorr.c
   @brief topological correlator measurement
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "momenta.h"
#include "resampled_ops.h"

//#define SUBZERO

// the the topological corr < Q >^2 / V
int
fit_Qcorr( struct input_params *Input )
{   
  size_t i , j , shift = 0 ;
  const double fac = 1.0/ ( 10*10*10*10. ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift = j ;
  }
  
  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}

// partially sum the correlator \sum_{x=0}^{r^2} < q(x+y) q(x) >
int
fit_Qsusc( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // subtract the zero by hand?
    #ifdef SUBZERO
    equate_constant( &Input -> Data.y[shift] , 0.0 ,
		     Input -> Data.y[shift].NSAMPLES ,
		     Input -> Data.y[shift].restype ) ;
    #endif
    
    // compute the partial sum up to whatever cutoff
    for( j = shift+1 ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      
      add( &Input -> Data.y[j] , Input -> Data.y[j-1] ) ;
    }
    
    shift = j ;
  }

  average_equivalent( Input ) ;

  // compute the derivative of the partial sum
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim; i++ ) {

    struct resampled der = init_dist( NULL ,
				      Input -> Data.x[shift].NSAMPLES ,
				      Input -> Data.x[shift].restype ) ;

    for( j = shift+1 ; j < shift + Input -> Data.Ndata[i] - 1 ; j++ ) {

	equate( &der , Input -> Data.y[j+1] ) ;
	subtract( &der , Input -> Data.y[j-1] ) ;
	mult_constant( &der , 1.0/( Input -> Data.x[j+1].avg - Input -> Data.x[j-1].avg ) ) ;
	fprintf( stdout , "[DER] %e %e %e \n" , Input -> Data.x[j].avg , der.avg , der.err ) ;

	if( ( der.avg - der.err ) < 0.0 ) {
	  fprintf( stdout , "CANDIDATE %e %e %e \n" , Input -> Data.x[j].avg ,
		   Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
	  break ;
	}
    }

    free( der.resampled ) ;
    
    shift = j ;
  }

  // multiply
  shift = 0 ;
  const double fac = 1.0 / ( 10*10*10*10. ) ;
  for( i = 0 ; i < Input -> Data.Nsim; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift = j ;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
	
  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}

// fit the slab method evaluation of Qsusc
int
fit_Qslab( struct input_params *Input )
{  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    const double fac = Input -> Traj[i].Dimensions[3] / ( 10 * 10 * 10 * 10. ) ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      add_constant( &Input -> Data.x[j] , 1.0 ) ;
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift=j;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
	
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}

// clean these up
#ifdef SUBZERO
  #undef SUBZERO
#endif
