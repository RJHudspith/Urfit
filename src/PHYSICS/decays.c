/**
   @file decays.c
   @brief compute decay constants from the fit
 */
#include "gens.h"

#include "resampled_ops.h"

// computes the decay constant for a given amplitude from our simultaneous fit
struct resampled
decay( const struct resampled *fitparams ,
       const struct input_params Input ,
       const size_t Mass_idx  ,
       const size_t Amp_idx )
{    
  struct resampled result = init_dist( &fitparams[ Amp_idx ] ,
				       fitparams[ Amp_idx ].NSAMPLES ,
				       fitparams[ Amp_idx ].restype ) ;

  const double vol_fac = 2./( Input.Traj[0].Dimensions[0] *
			      Input.Traj[0].Dimensions[1] *
			      Input.Traj[0].Dimensions[2] ) ;

  // f = sqrt{ 2/V*(A^L)^2/m_\pi }
  mult( &result , fitparams[ Amp_idx ] ) ;
  mult_constant( &result , vol_fac ) ;
  divide( &result , fitparams[ Mass_idx ] ) ;
  root( &result ) ;

  fprintf( stdout , "Decay/Z_A :: %e,%e \n" , result.avg , result.err ) ;
  
  return result ;
}
