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

static double
get_t0MK( const double aM )
{
  return 1./(1.8767 + 1.8333*aM + 0.12686*aM*aM) ;
}

int
su2_shit( struct input_params *Input )
{  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    equate_constant( &Input -> Data.x[i] ,
		     get_t0MK( Input->Data.x[i].avg ) ,
		     Input -> Data.x[i].NSAMPLES ,
		     Input -> Data.x[i].restype ) ;
    mult_constant( &Input -> Data.y[i] , 1.3317 ) ;
  }

  // finite volume
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift+Input->Data.Ndata[i] ; j++ ) {
      fprintf( stdout , "%e %1.15e %1.15e\n" ,
	       1.3317/(Input -> Traj[i].Dimensions[0]) ,
	       Input -> Data.y[j].avg ,
	       Input -> Data.y[j].err ) ;
    }
    shift += Input->Data.Ndata[i] ;
    printf( "\n" ) ;
  }

  double chisq ;  
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
