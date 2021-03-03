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
#include "plot_fitfunc.h"
#include "stats.h"

static double
get_t0MK( const double aM )
{
  return 1./(1.8767 + 1.8333*aM + 0.12686*aM*aM) ;
}

int
su2_shit( struct input_params *Input )
{
  /*
  raise( &Input -> Data.y[0] , 0.25 ) ;
  printf( "U0 %f %f\n" , Input ->Data.y[0].avg , Input -> Data.y[0].err ) ;
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
  */

  /*
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift+Input->Data.Ndata[i] ; j++ ) {
      //raise( &Input -> Data.y[j] , 2 ) ;
      fprintf( stdout , "%e %1.15e %1.15e\n" ,
	       Input -> Data.x[j].avg ,
	       Input -> Data.y[j].avg ,
	       Input -> Data.y[j].err ) ;
    }
	shift = j ;
  }
  */
  const size_t end = Input -> Data.Ndata[0]-1 ;
  
  struct resampled sub = init_dist( &Input -> Data.x[ end ] , Input -> Data.x[ end ].NSAMPLES , Input -> Data.x[ end ].restype ) ;
  for( int i = 0 ; i < Input -> Data.Ndata[0] ; i++ ) {
    //root( &Input -> Data.y[i] ) ;
    subtract( &Input -> Data.x[i] , sub ) ;
  }
  free( sub.resampled ) ;
  
  double chisq ;  
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
