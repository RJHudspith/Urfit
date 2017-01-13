/**
   @file fit_chooser.c
   @brief choose and set up the fit
 */
#include "gens.h"

#include "fits.h"      // is a list of fit headers
#include "ffunction.h"
#include "pmap.h"      // for allocating the pmap

// return number of params in the fit
size_t
get_Nparam( const struct fit_info Fit )
{
 switch( Fit.Fitdef ) {
 case COSH : return 2 ;
 case EXP : return 2 ;
 case EXP_PLUSC : return 3 ;
 case PADE : return Fit.N + Fit.M + 1 ;
 case POLY : return Fit.N + 1 ;
 case SINH : return 2 ;
 }
 return 0 ;
}

// initialise the fit
struct fit_descriptor
init_fit( const struct data_info Data ,
	  const struct fit_info Fit )
{
  struct fit_descriptor fdesc ;
  switch( Fit.Fitdef ) {
  case COSH : 
    fdesc.func       = fcosh ;
    fdesc.F          = cosh_f ;
    fdesc.dF         = cosh_df ;
    fdesc.d2F        = cosh_d2f ;
    fdesc.guesses    = cosh_guesses ;
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
    fdesc.guesses    = exp_plusc_guesses ;
    break ;
  case PADE :
    pade_set_nm( Fit.N , Fit.M ) ;
    fdesc.func       = fpade ;
    fdesc.F          = pade_f ;
    fdesc.dF         = pade_df ;
    fdesc.d2F        = pade_d2f ;
    fdesc.guesses    = pade_guesses ;
    break ;
  case POLY :
    poly_set_n( Fit.N + 1 ) ;
    fdesc.func       = fpoly ;
    fdesc.F          = poly_f ;
    fdesc.dF         = poly_df ;
    fdesc.d2F        = poly_d2f ;
    fdesc.guesses    = poly_guesses ;
    break ;
  case SINH : 
    fdesc.func       = fsinh ;
    fdesc.F          = sinh_f ;
    fdesc.dF         = sinh_df ;
    fdesc.d2F        = sinh_d2f ;
    fdesc.guesses    = sinh_guesses ;
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
    //Data.Nsim * ( fdesc.Nparam - Ncommon ) + Ncommon ;

  // allocate the fitfunction
  fdesc.f = allocate_ffunction( fdesc.Nlogic , Data.Ntot ) ;
  
  // corrfit
  fdesc.f.CORRFIT = Fit.Corrfit ;

  return fdesc ;
}
