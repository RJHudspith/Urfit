/**
   @file correlator.c
   @brief correlator analysis
 */
#include "gens.h"

#include "decays.h"
#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"

int
correlator_analysis( struct input_params *Input )
{
  // compute an effective mass ? TODO
  struct resampled *effmass = effective_mass( Input , ASINH_EFFMASS ) ;

  size_t i ;
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  // compute a decay constant
  if( Input -> Fit.Fitdef == PP_AA_WW ||
      Input -> Fit.Fitdef == PP_AA ) {
    decay( Fit , *Input ) ;
  }

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
