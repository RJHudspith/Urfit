/**
   @file HLBL.c
   @brief integrates the HLBL from Jamie's data
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "Nint.h"
#include "resampled_ops.h"

int
HLBL_analysis( struct input_params *Input )
{
  const double ammap[7] = { 0.375 , 0.25 , 0.1875 , 0.125 , 0.09375 , 0.083333 , 0.0625 } ;
  
  const double alpha_QED = 1.0/137.035999 ;
  const double prefac = pow(4*M_PI*alpha_QED,3)/3 ;
  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // multiply by |x|^3*prefac
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      // x is y^3 when comparing to note
      root( &Input -> Data.x[j] ) ;

      mult_constant( &Input -> Data.y[j] ,
		     ammap[i]*prefac*pow( Input -> Data.x[j].avg , 3 ) ) ;
    }

    fprintf( stdout , "\nPerforming integrals\n" ) ; 
    size_t k ;
    for( k = 2 ; k < Input -> Data.Ndata[i] ; k++ ) {
      // do the numberical integration
      struct resampled Int = Nint( Input -> Data.x + shift ,
				   Input -> Data.y + shift ,
				   k , true ) ;
      fprintf( stdout , "Integral %zu %e %e %e \n" ,
	       i , Input -> Data.x[shift+k].avg , Int.avg , Int.err ) ;
      free( Int.resampled ) ;
    }

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      //mult_constant( &Input -> Data.x[j] , ammap[i] ) ;
    }
    
    shift += Input -> Data.Ndata[i] ;
  }


  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot_and_Nint( *Input , &Chi ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;  
  
  return SUCCESS ;
}
