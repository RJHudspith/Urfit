/**
   @file fit_chooser.c
   @brief choose and set up the fit
 */
#include "gens.h"

#include "cosh.h"
#include "exp.h"
#include "exp_plusc.h"
#include "pade.h"
#include "poly.h"
#include "pp_aa.h"
#include "pp_aa_ww.h"
#include "sinh.h"

#include "ffunction.h"
#include "pmap.h"      // for allocating the pmap

// return number of params in the fit
size_t
get_Nparam( const struct fit_info Fit )
{
 switch( Fit.Fitdef ) {
   // fall through as they are all the same, multi-exp,cosh or sinh
 case COSH :
 case SINH :
 case EXP :
   return 2 * Fit.N ;
   // these have independent amounts
 case EXP_PLUSC : return 3 ;
 case PADE : return Fit.N + Fit.M ;
 case POLY : return Fit.N + 1 ;
 case PP_AA : return 5 ;
 case PP_AA_WW : return 5 ;
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
  case SINH : 
    fdesc.func       = fsinh ;
    fdesc.F          = sinh_f ;
    fdesc.dF         = sinh_df ;
    fdesc.d2F        = sinh_d2f ;
    fdesc.guesses    = exp_guesses ;
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
