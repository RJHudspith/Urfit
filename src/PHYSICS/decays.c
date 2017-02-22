/**
   @file decays.c
   @brief compute decay constants from the fit
 */
#include "gens.h"

#include "resampled_ops.h"

// computes the decay constant for a given amplitude from our simultaneous fit
struct resampled
decay( const struct resampled *fitparams ,
       const struct input_params Input )
{
  enum { MASS = 0 , PL = 1 , AL = 2 , PW = 3 , AW = 4 }  ;
    
  struct resampled result ;
  result.resampled = malloc( fitparams[0].NSAMPLES * sizeof( double ) ) ;

  const double vol_fac = 2.0 / ( Input.Traj[0].Dimensions[0] *
				 Input.Traj[0].Dimensions[1] *
				 Input.Traj[0].Dimensions[2] ) ;

  // f = 
  equate( &result , fitparams[ PL ] ) ;
  mult( &result , fitparams[ PL ] ) ;
  mult_constant( &result , vol_fac ) ;
  divide( &result , fitparams[ MASS ] ) ;
  root( &result ) ;

  printf( "Decay :: %e +/- %e \n" , result.avg , result.err ) ;
  
  return result ;
}
