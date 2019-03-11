/**
   @file su2_shit.c
   @brief as the file name suggests I wasn't fond of this analysis
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "fit_chooser.h"
#include "ffunction.h"
#include "init.h"
#include "make_xmgrace.h"
#include "resampled_ops.h"
#include "stats.h"

int
su2_shit( struct input_params *Input )
{
  size_t i = 0 ;
  double chisq ;

  const double xmap[10] = { 0.2606 , 0.2036 , 0.2290 ,
			    0.1827 , 0.1657 , 0.1179 ,
			    0.0793 , 0.0578 , 0.0444 ,
			    0.0353 } ;
  
  size_t j ;
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    equate_constant( &Input -> Data.x[i] , xmap[i%10] ,
		     Input -> Data.x[i].NSAMPLES ,
		     Input -> Data.x[i].restype ) ;
    mult_constant( &Input -> Data.y[i] , 1.3317 ) ;
  }

  for( j = 0 ; j < 10 ; j++ ) {
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      fprintf( stdout , "%e %1.15e %1.15e\n" ,
	       1.3317*1.3317*1.3317/(Input -> Traj[i].Dimensions[0]*
		   Input -> Traj[i].Dimensions[1]*
		   Input -> Traj[i].Dimensions[2]) ,
	       Input -> Data.y[j+10*i].avg ,
	       Input -> Data.y[j+10*i].err ) ;
    }
    printf( "\n" ) ;
  }
  
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
