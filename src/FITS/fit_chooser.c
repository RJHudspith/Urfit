/**
   @file fit_chooser.c
   @brief choose and set up the fit
 */
#include "gens.h"

#include "alpha_D0.h"
#include "alpha_D0_multi.h"
#include "cosh.h"
#include "exp.h"
#include "exp_plusc.h"
#include "pade.h"
#include "poly.h"
#include "pp_aa.h"
#include "pp_aa_ww.h"
#include "ppaa.h"
#include "poles.h"
#include "Qcorr_bessel.h"
#include "qslab.h"
#include "sinh.h"

#include "ffunction.h"
#include "pmap.h"      // for allocating the pmap

// return number of params in the fit
size_t
get_Nparam( const struct fit_info Fit )
{
 switch( Fit.Fitdef ) {
 case ALPHA_D0 :
   return 4 ;
 case ALPHA_D0_MULTI :
   return 4 ;
   // fall through as they are all the same, multi-exp,cosh or sinh
 case COSH :
 case SINH :
 case EXP :
   return 2 * Fit.N ;
   // these have independent amounts
 case EXP_PLUSC : return 3 ;
 case PADE : return Fit.N + Fit.M ;
 case POLES : return 7 ;
 case POLY : return Fit.N + 1 ;
 case PPAA : return 3 ;
 case PP_AA : return 5 ;
 case PP_AA_WW : return 5 ;
 case QCORR_BESSEL : return 2 ;
 case QSLAB : return 3 ;
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
    fdesc.func       = falpha_D0_multi ;
    fdesc.F          = alpha_D0_multi_f ;
    fdesc.dF         = alpha_D0_multi_df ;
    fdesc.d2F        = alpha_D0_multi_d2f ;
    fdesc.guesses    = alpha_D0_multi_guesses ;
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
    //fdesc.guesses    = exp_guesses ;
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
  case PP_AA_WW :
    fdesc.func       = fpp_aa_ww ;
    fdesc.F          = pp_aa_ww_f ;
    fdesc.dF         = pp_aa_ww_df ;
    fdesc.d2F        = pp_aa_ww_d2f ;
    fdesc.guesses    = pp_aa_ww_guesses ;
    break ;
  case POLES :
    fdesc.func       = fpoles ;
    fdesc.F          = poles_f ;
    fdesc.dF         = poles_df ;
    fdesc.d2F        = poles_d2f ;
    fdesc.guesses    = poles_guesses ;
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

  // allocate the fitfunction
  fdesc.f = allocate_ffunction( fdesc.Nlogic , Data.Ntot ) ;
  
  // corrfit
  fdesc.f.CORRFIT = Fit.Corrfit ;

  return fdesc ;
}
