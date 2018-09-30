/**
   @file fit_chooser.c
   @brief choose and set up the fit
 */
#include "gens.h"

#include "alpha_D0.h"
#include "alpha_D0_multi.h"
#include "alpha_D0_multi_new.h"
#include "adler_alpha_D0.h"
#include "adler_alpha_D0_multi.h"
#include "cornell.h"
#include "cosh.h"
#include "exp.h"
#include "exp_plusc.h"
#include "fvol1.h"
#include "c4c7.h"
#include "nrqcd_exp.h"
#include "nrqcd_exp2.h"
#include "pade.h"
#include "poly.h"
#include "pp_aa.h"
#include "pp_aa_exp.h"
#include "pp_aa_ww.h"
#include "pp_aa_ww_r2.h"
#include "ppaa.h"
#include "poles.h"
#include "Qcorr_bessel.h"
#include "qsusc_su2.h"
#include "qslab.h"
#include "sinh.h"
#include "udcb_heavy.h"

#include "ffunction.h"
#include "pmap.h"      // for allocating the pmap
#include "su2_shitfit.h"

// return number of params in the fit
size_t
get_Nparam( const struct fit_info Fit )
{
 switch( Fit.Fitdef ) {
 case ALPHA_D0 :
   return 4 ;
 case ALPHA_D0_MULTI :
   //return 4 ;
   return 6 ;
 case ADLERALPHA_D0 :
   return 4 ;
 case ADLERALPHA_D0_MULTI :
   return 4 ;
   // fall through as they are all the same, multi-exp,cosh or sinh
 case C4C7 :
   return 3 ;
 case CORNELL : return 6 ;
 case COSH :
 case SINH :
 case EXP :
   return 2 * Fit.N ;
   // these have independent amounts
 case EXP_PLUSC : return 3 ;
 case FVOL1 : return 2 ;
 case PADE : return Fit.N + Fit.M ;
 case POLES : return Fit.N + Fit.M + 1 ;
 case POLY : return Fit.N + 1 ;
 case PPAA : return 3*Fit.N ;
 case PP_AA : return 5 ;
 case PP_AA_EXP : return 5 ;
 case PP_AA_WW : return 5 ;
 case PP_AA_WW_R2 : return 5 + 2*Fit.N ;
 case QCORR_BESSEL : return 2 ;
 case QSLAB : return 3 ;
 case QSUSC_SU2 : return 3 ;
 case UDCB_HEAVY : return 5 ;
 case NRQCD_EXP : return 2 ;
 case NRQCD_EXP2 : return 5 ;
 case SU2_SHITFIT : return 2 ;
 case NOFIT : return 0 ;
 }
 return 0 ;
}

// initialise the fit
struct fit_descriptor
init_fit( const struct data_info Data ,
	  const struct fit_info Fit )
{
  struct fit_descriptor fdesc ;
  fdesc.linmat = NULL ;
  
  switch( Fit.Fitdef ) {
  case ALPHA_D0 :
    fdesc.func       = falpha_D0 ;
    fdesc.F          = alpha_D0_f ;
    fdesc.dF         = alpha_D0_df ;
    fdesc.d2F        = alpha_D0_d2f ;
    fdesc.guesses    = alpha_D0_guesses ;
    break ;
  case ALPHA_D0_MULTI :
    /*
    fdesc.func       = falpha_D0_multi ;
    fdesc.F          = alpha_D0_multi_f ;
    fdesc.dF         = alpha_D0_multi_df ;
    fdesc.d2F        = alpha_D0_multi_d2f ;
    fdesc.guesses    = alpha_D0_multi_guesses ;
    */
    fdesc.func       = falpha_D0_multi2 ;
    fdesc.F          = alpha_D0_multi2_f ;
    fdesc.dF         = alpha_D0_multi2_df ;
    fdesc.d2F        = alpha_D0_multi2_d2f ;
    fdesc.guesses    = alpha_D0_multi2_guesses ;
    break ;
  case ADLERALPHA_D0 :
    fdesc.func       = fadleralpha_D0 ;
    fdesc.F          = adleralpha_D0_f ;
    fdesc.dF         = adleralpha_D0_df ;
    fdesc.d2F        = adleralpha_D0_d2f ;
    fdesc.guesses    = adleralpha_D0_guesses ;
    break ;
  case ADLERALPHA_D0_MULTI :
    fdesc.func       = fadleralpha_D0_multi ;
    fdesc.F          = adleralpha_D0_multi_f ;
    fdesc.dF         = adleralpha_D0_multi_df ;
    fdesc.d2F        = adleralpha_D0_multi_d2f ;
    fdesc.guesses    = adleralpha_D0_multi_guesses ;
    break ;
  case C4C7 :
    fdesc.func       = fc4c7 ;
    fdesc.F          = c4c7_f ;
    fdesc.dF         = c4c7_df ;
    fdesc.d2F        = c4c7_d2f ;
    fdesc.guesses    = c4c7_guesses ;
    break ;
  case CORNELL :
    fdesc.func       = fcornell ;
    fdesc.F          = cornell_f ;
    fdesc.dF         = cornell_df ;
    fdesc.d2F        = cornell_d2f ;
    fdesc.guesses    = cornell_guesses ;
    break ;
  case COSH : 
    fdesc.func       = fcosh ;
    fdesc.F          = cosh_f ;
    fdesc.dF         = cosh_df ;
    fdesc.d2F        = cosh_d2f ;
    fdesc.guesses    = exp_guesses ;
    break ;
  case EXP :
    fdesc.func       = fexp ;
    fdesc.F          = exp_f ;
    fdesc.dF         = exp_df ;
    fdesc.d2F        = exp_d2f ;
    fdesc.guesses    = exp_guesses ; 
    break ;
  case EXP_PLUSC : 
    fdesc.func       = fexp_plusc ;
    fdesc.F          = exp_plusc_f ;
    fdesc.dF         = exp_plusc_df ;
    fdesc.d2F        = exp_plusc_d2f ;
    break ;
  case FVOL1 :
    fdesc.func       = ffvol1 ;
    fdesc.F          = fvol1_f ;
    fdesc.dF         = fvol1_df ;
    fdesc.d2F        = fvol1_d2f ;
    break ;
  case NRQCD_EXP :
    fdesc.func       = fnrqcd_exp ;
    fdesc.F          = nrqcd_exp_f ;
    fdesc.dF         = nrqcd_exp_df ;
    fdesc.d2F        = nrqcd_exp_d2f ;
    fdesc.guesses    = nrqcd_exp_guesses ; 
    break ;
  case NRQCD_EXP2 :
    fdesc.func       = fnrqcd_exp2 ;
    fdesc.F          = nrqcd_exp2_f ;
    fdesc.dF         = nrqcd_exp2_df ;
    fdesc.d2F        = nrqcd_exp2_d2f ;
    fdesc.guesses    = nrqcd_exp2_guesses ; 
    break ;
  case PADE :
    fdesc.func       = fpade ;
    fdesc.F          = pade_f ;
    fdesc.dF         = pade_df ;
    fdesc.d2F        = pade_d2f ;
    fdesc.guesses    = pade_guesses ;
    break ;
  case POLY :
    fdesc.func       = fpoly ;
    fdesc.F          = poly_f ;
    fdesc.dF         = poly_df ;
    fdesc.d2F        = poly_d2f ;
    fdesc.guesses    = poly_guesses ;
    fdesc.linmat     = poly_linmat ;
    break ;
  case PPAA :
    fdesc.func       = fppaa ;
    fdesc.F          = ppaa_f ;
    fdesc.dF         = ppaa_df ;
    fdesc.d2F        = ppaa_d2f ;
    fdesc.guesses    = ppaa_guesses ;
    break ;
  case PP_AA :
    fdesc.func       = fpp_aa ;
    fdesc.F          = pp_aa_f ;
    fdesc.dF         = pp_aa_df ;
    fdesc.d2F        = pp_aa_d2f ;
    fdesc.guesses    = pp_aa_guesses ;
    break ;
  case PP_AA_EXP :
    fdesc.func       = fpp_aa_exp ;
    fdesc.F          = pp_aa_exp_f ;
    fdesc.dF         = pp_aa_exp_df ;
    fdesc.d2F        = pp_aa_exp_d2f ;
    fdesc.guesses    = pp_aa_exp_guesses ;
    break ;
  case PP_AA_WW :
    fdesc.func       = fpp_aa_ww ;
    fdesc.F          = pp_aa_ww_f ;
    fdesc.dF         = pp_aa_ww_df ;
    fdesc.d2F        = pp_aa_ww_d2f ;
    fdesc.guesses    = pp_aa_ww_guesses ;
    break ;
  case PP_AA_WW_R2 :
    fdesc.func       = fpp_aa_ww_r2 ;
    fdesc.F          = pp_aa_ww_r2_f ;
    fdesc.dF         = pp_aa_ww_r2_df ;
    fdesc.d2F        = pp_aa_ww_r2_d2f ;
    fdesc.guesses    = pp_aa_ww_r2_guesses ;
    break ;
  case POLES :
    fdesc.func       = fpoles ;
    fdesc.F          = poles_f ;
    fdesc.dF         = poles_df ;
    fdesc.d2F        = poles_d2f ;
    fdesc.guesses    = poles_guesses ;
    fdesc.linmat     = poles_linmat ;
    break ;
  case SINH : 
    fdesc.func       = fsinh ;
    fdesc.F          = sinh_f ;
    fdesc.dF         = sinh_df ;
    fdesc.d2F        = sinh_d2f ;
    fdesc.guesses    = exp_guesses ;
    break ;
  case QCORR_BESSEL :
    fdesc.func       = fQcorr_bessel ;
    fdesc.F          = Qcorr_bessel_f ;
    fdesc.dF         = Qcorr_bessel_df ;
    fdesc.d2F        = Qcorr_bessel_d2f ;
    fdesc.guesses    = Qcorr_bessel_guesses ;
    break ;
  case QSLAB :
    fdesc.func       = fqslab ;
    fdesc.F          = qslab_f ;
    fdesc.dF         = qslab_df ;
    fdesc.d2F        = qslab_d2f ;
    fdesc.guesses    = qslab_guesses ;
    break ;
  case QSUSC_SU2 :
    fdesc.func       = fqsusc_su2 ;
    fdesc.F          = qsusc_su2_f ;
    fdesc.dF         = qsusc_su2_df ;
    fdesc.d2F        = qsusc_su2_d2f ;
    fdesc.guesses    = qsusc_su2_guesses ;
    break ;
  case UDCB_HEAVY :
    fdesc.func       = fudcb_heavy ;
    fdesc.F          = udcb_heavy_f ;
    fdesc.dF         = udcb_heavy_df ;
    fdesc.d2F        = udcb_heavy_d2f ;
    fdesc.guesses    = udcb_heavy_guesses ;
    break ;
  case SU2_SHITFIT :
    fdesc.func       = fsu2_shitfit ;
    fdesc.F          = su2_shitfit_f ;
    fdesc.dF         = su2_shitfit_df ;
    fdesc.d2F        = su2_shitfit_d2f ;
    fdesc.guesses    = su2_shitfit_guesses ;
    break ;
  case NOFIT :
    break ;
  }

  fdesc.Nparam = get_Nparam( Fit ) ;

  // count how many common parameters we have
  size_t Ncommon = 0 , i ;
  for( i = 0 ; i < fdesc.Nparam ; i++ ) {
    if( Fit.Sims[i] == true ) {
      Ncommon++ ;
    }
  }
  fdesc.Nlogic = Fit.Nlogic ;

  // set the N and Ms of the fit
  fdesc.N = Fit.N ;
  fdesc.M = Fit.M ;

  // allocate the fitfunction
  fdesc.f = allocate_ffunction( fdesc.Nlogic , Data.Ntot ) ;
  
  // corrfit
  fdesc.f.CORRFIT = Fit.Corrfit ;

  return fdesc ;
}
