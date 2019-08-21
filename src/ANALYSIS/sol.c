/**
   @brief compute the speed of light from

   (E(p)^2-E(0)^2)/p^2 = c^2
 */
#include "gens.h"

#include "resampled_ops.h"

static double
psq( struct input_params *Input , const size_t i )
{
  double psq = 0.0 ;
  size_t mu ;
  for( mu = 0 ; mu < 3 ; mu++ ) {
    double x = Input -> Traj[i].mom[mu]*2.*M_PI/Input->Traj[i].Dimensions[mu] ;
    psq += x*x ;
  }
  return psq ;
}

int
sol_analysis( struct input_params *Input )
{
  printf( "Here\n" ) ;
  size_t i ;
  raise( &Input -> Data.y[0] , 2 ) ;
  for( i = Input -> Data.Nsim-1 ; i > 0 ; i-- ) {
    raise( &Input -> Data.y[i] , 2 ) ;
    subtract( &Input -> Data.y[i] , Input -> Data.y[0] ) ;
    const double q2 = psq( Input , i ) ;
    divide_constant(  &Input -> Data.y[i] ,  q2 ) ;
    printf( "%f %f %f\n" , q2 , Input -> Data.y[i].avg , Input -> Data.y[i].err ) ;
  }
  
  return SUCCESS ;
}
