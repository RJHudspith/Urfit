#include "gens.h"

#include "decays.h"
#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"

//#define FIT_EFFMASS
//#define COMBINE

int
wall_local_analysis( struct input_params *Input )
{
  if( Input -> Data.Nsim%2 == 1 ) {
    fprintf( stderr , "[WL] must have an even number of correlators\n" ) ;
    return FAILURE ;
  }
  
  const size_t NWL = Input -> Data.Nsim / 2 ;
  const size_t LT = Input -> Traj[0].Dimensions[3] / 2 ;

  size_t i , j ;
  for( i = 0 ; i < NWL ; i++ ) {
    for( j = LT*i ; j < (i+1)*LT ; j++ ) {
      raise( &(Input -> Data.y[j]) , 2 ) ;
      divide( &(Input -> Data.y[j]) ,
	      Input -> Data.y[j+NWL*LT] ) ;
    }
  }

  printf( "Effmass\n" ) ;
    
  // compute an effective mass ? TODO
  struct resampled *effmass = effective_mass( Input , ASINH_EFFMASS ) ;

#ifdef FIT_EFFMASS
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    equate( &Input -> Data.y[i] , effmass[i] ) ;
  }
#endif

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}
