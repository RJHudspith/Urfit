/**
   @brief quick check of the exceptional fits
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "resampled_ops.h"

#include "init.h"

#define DOUBLE_POLE
//#define SINGLE_POLE
//#define LINEAR_FIT

int
fit_exceptional( struct input_params *Input )
{
  #ifdef LINEAR_FIT
  if( Input -> Fit.Fitdef != POLY && Input -> Fit.N == 1 ) {
    fprintf( stderr , "Linear fit expects linear polynomial fit\n" ) ;
    return FAILURE ;
  }
  #endif
  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      #ifdef DOUBLE_POLE
      mult_constant( &Input -> Data.y[j] ,
		     Input -> Data.x[j].avg * Input -> Data.x[j].avg ) ;
      #elif ( defined SINGLE_POLE ) || (defined  LINEAR_FIT)
      mult_constant( &Input -> Data.y[j] , Input -> Data.x[j].avg ) ;
      #endif
      if( i == 1 ) {
	mult_constant( &Input -> Data.y[j] , -1 ) ;
      }
    }
    shift = j ;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  if( Input -> Fit.Fitdef == PADE &&
      Input -> Fit.N == 2 &&
      Input -> Fit.M == 1 ) {
    fprintf( stdout , "[EX] In here \n" ) ;
      
    struct resampled temp = init_dist( &fit[0] , fit[0].NSAMPLES ,
				       fit[0].restype ) ;

    fprintf( stdout , "[EX] Coeff 0 : %f %f \n" , temp.avg , temp.err ) ;

    // c1 = \pi_1 - \pi_0 * \pi_2
    mult( &temp , fit[2] ) ;
    subtract( &temp , fit[1] ) ;
    mult_constant( &temp , -1 ) ;

    fprintf( stdout , "[EX] Coeff 1 : %f %f \n" , temp.avg , temp.err ) ;
      
    // c2 = -c1 * pi_2
    mult( &temp , fit[2] ) ;
    mult_constant( &temp , -1 ) ;

    fprintf( stdout , "[EX] Coeff 2 : %f %f \n" , temp.avg , temp.err ) ;

    // c3 = -c2 * pi_2
    mult( &temp , fit[2] ) ;
    mult_constant( &temp , -1 ) ;
      
    fprintf( stdout , "[EX] Coeff 3 : %f %f \n" , temp.avg , temp.err ) ;
      
    free( temp.resampled ) ;
  }

  #ifdef LINEAR_FIT
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    // subtract the pole
    subtract( &Input -> Data.y[i] , fit[0] ) ;
    // divide data through by m again
    mult_constant( &Input -> Data.y[i] , 1.0 / ( Input -> Data.x[i].avg ) ) ;
  }
  #endif

  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  #ifdef LINEAR_FIT
  fit = fit_and_plot( *Input , &chisq ) ;
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  #endif
  
  return SUCCESS ;
}
