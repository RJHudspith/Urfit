/**
   @file an_wrapper.c
   @brief analysis function wrapper
 */
#include "gens.h"

#include "alphas.h"
#include "beta_crit.h"
#include "correlator.h"
#include "fit_and_plot.h"
#include "exceptional.h"
#include "hvp_pade.h"
#include "Qcorr.h"
#include "sun_flow.h"

int
an_wrapper( struct input_params *Input )
{
  double chi = 0.0 ;
  switch( Input -> Analysis ) {
  case Alphas :
    return fit_alphas( Input ) ;
  case Beta_crit :
    return beta_crit( Input ) ;
  case Correlator :
    return correlator_analysis( Input ) ;
  case Exceptional :
    return fit_exceptional( Input ) ;
  case HVP :
    return fit_hvp( Input ) ;
  case Qcorr :
    return fit_Qcorr( Input ) ;
  case Qsusc :
    return fit_Qsusc( Input ) ;
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
