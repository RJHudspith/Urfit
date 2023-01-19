/**
   @file an_wrapper.c
   @brief analysis function wrapper
 */
#include "gens.h"

#include "adler.h"
#include "alphas.h"
#include "beta_crit.h"
#include "binding_corr.h"
#include "correlator.h"
#include "CPCV.h"
#include "fit_and_plot.h"
#include "exceptional.h"
#include "general_ops.h"
#include "HLBL.h"
#include "hvp_pade.h"
#include "inverse.h"
#include "KKops.h"
#include "KK_BK.h"
#include "nrqcd.h"
#include "nrqcd_old.h"
#include "nrqcd_slope.h"
#include "Qcorr.h"
#include "Qcorr_fixed.h"
#include "Qmoments.h"
#include "Ren_rats.h"
#include "statpot.h"
#include "sun_flow.h"
#include "tetra_gevp.h"
#include "Wall_Local.h"
#include "ZV.h"

#include "nrqcd_baremass.h"
#include "c4c7_analysis.h"
#include "su2_shit.h"
#include "pof.h"
#include "sol.h"
#include "SpinOrbit.h"

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
  case Binding_Corr :
    return binding_corr_analysis2( Input ) ;
  case Correlator :
    return correlator_analysis( Input ) ;
  case Exceptional :
    return fit_inverse( Input ) ;
  case General :
    return gen_ops( Input ) ;
  case HLBL :
    return HLBL_analysis( Input ) ;
  case HVP :
    return fit_hvp( Input ) ;
  case KKops :
    return fit_ratios( Input ) ;
  case KK_BK :
    return fit_bags( Input ) ;
  case Nrqcd :
    return nrqcd_analysis( Input ) ;
    //return nrqcd_old_analysis( Input ) ;
    //return nrqcd_slope_analysis( Input ) ;
  case Pof :
    return pof_analysis( Input ) ;
    //return pof_analysis_fixed( Input ) ;
  case Qcorr :
    //return fit_Qcorr( Input ) ;
    return Qmoments( Input ) ;
    //return TraditionalQ( Input ) ;
    //return CumFromMom( Input ) ;
  case Qsusc :
    return fit_Qsusc( Input ) ;
  case Qslab :
    return fit_Qslab( Input ) ;
  case QslabFix :
    return fit_Qslab_fixed( Input ) ;
  case Ren_Rats :
    return renormalise_rats( Input ) ;
  case SpinOrbit :
    return spin_orbit( Input ) ;
  case Sol :
    return sol_analysis( Input ) ;
  case StaticPotential :
    return statpot_analysis_v2( Input ) ;
  case TetraGEVP :
    return tetra_gevp_analysis( Input ) ;
  case TetraGEVP_Fixed :
    return tetra_gevp_fixed_delta_analysis( Input ) ;
  case Wflow :
    //return sun_wflow_analysis( Input ) ;
    return sun_set( Input ) ;
  case ZV :
    return ZV_analysis( Input ) ;
  case Fit :
    return nrqcd_baremass_analysis( Input ) ;
    //return c4c7_analysis( Input ) ;
    //return su2_shit( Input ) ;
    //return HAL_analysis( Input ) ;
    //return sol_analysis( Input ) ;
    //return binding_corr_analysis2( Input ) ;
  default :
    if( fit_and_plot( *Input , &chi ) == NULL ) {
      return FAILURE ;
    }
    break ;
  }
 return SUCCESS ;
}
