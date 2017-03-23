/**
   @file an_wrapper.c
   @brief analysis function wrapper
 */
#include "gens.h"

#include "alphas.h"
#include "correlator.h"
#include "fit_and_plot.h"
#include "hvp_pade.h"
#include "sun_flow.h"

int
an_wrapper( struct input_params *Input )
{
  double chi = 0.0 ;
  switch( Input -> Analysis ) {
  case Alphas :
    return fit_alphas( Input ) ;
  case Correlator :
    return correlator_analysis( Input ) ;
  case HVP :
    return fit_hvp( Input ) ;
  case Wflow :
    return sun_wflow_analysis( Input ) ;
  case Fit :
  default :
    if( fit_and_plot( *Input , &chi ) == NULL ) {
      return FAILURE ;
    }
    break ;
  }
 return SUCCESS ;
}
