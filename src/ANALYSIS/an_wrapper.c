/**
   @file an_wrapper.c
   @brief analysis function wrapper
 */
#include "gens.h"

#include "adler.h"
#include "alphas.h"
#include "beta_crit.h"
#include "correlator.h"
#include "CPCV.h"
#include "fit_and_plot.h"
#include "exceptional.h"
#include "general_ops.h"
#include "hvp_pade.h"
#include "inverse.h"
#include "KKops.h"
#include "KK_BK.h"
#include "Qcorr.h"
#include "Qmoments.h"
#include "Ren_rats.h"
#include "statpot.h"
#include "sun_flow.h"
#include "tetra_gevp.h"
#include "Wall_Local.h"

int
an_wrapper( struct input_params *Input )
{
  double chi = 0.0 ;
  switch( Input -> Analysis ) {
  case Adler :
    return adler_analysis( Input ) ;
  case Alphas :
    return fit_alphas( Input ) ;
  case Beta_crit :
    return beta_crit( Input ) ;
  case Correlator :
    return correlator_analysis( Input ) ;
    //return CVCP_analysis( Input ) ;
    //return wall_local_analysis( Input ) ;
  case Exceptional :
    //return fit_exceptional( Input ) ;
    return fit_inverse( Input ) ;
  case General :
    return gen_ops( Input ) ;
  case HVP :
    return fit_hvp( Input ) ;
  case KKops :
    return fit_ratios( Input ) ;
  case KK_BK :
    return fit_bags( Input ) ;
  case Qcorr :
    //return fit_Qcorr( Input ) ;
    return Qmoments( Input ) ;
  case Qsusc :
    return fit_Qsusc( Input ) ;
  case Qslab :
    return fit_Qslab( Input ) ;
  case Ren_Rats :
    return renormalise_rats( Input ) ;
  case StaticPotential :
    return statpot_analysis( Input ) ;
  case TetraGEVP :
    return tetra_gevp_analysis( Input ) ;
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
