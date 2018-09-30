#include "gens.h"

#include "fit_and_plot.h"

// just a linear fit
int
nrqcd_slope_analysis( struct input_params *Input )
{
  size_t i ;
  
  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  // invert the slope
  raise( &fit[1] , -1. ) ;
  mult_constant( &fit[1] , 0.5 ) ;
  
  printf( "MKIN %e %e \n" , fit[1].avg , fit[1].err ) ;
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
