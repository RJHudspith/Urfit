/**
   @file fit_chooser.c
   @brief choose and set up the fit
 */
#include "gens.h"

#include "fits.h" // is a list of fit headers
#include "ffunction.h"

// initialise the fit
struct fit_descriptor
init_fit( const fittype fit ,
	  const size_t Ndata ,
	  const corrtype CORRFIT )
{
  struct fit_descriptor fdesc ;
  switch( fit ) {
  case COSH : 
    fdesc.func       = fcosh ;
    fdesc.F          = cosh_f ;
    fdesc.dF         = cosh_df ;
    fdesc.d2F        = cosh_d2f ;
    fdesc.guesses    = cosh_guesses ; 
    fdesc.set_priors = cosh_priors ;
    fdesc.NPARAMS    = 2 ;
    break ;
  case EXP : 
    fdesc.func       = fexp ;
    fdesc.F          = exp_f ;
    fdesc.dF         = exp_df ;
    fdesc.d2F        = exp_d2f ;
    fdesc.guesses    = exp_guesses ; 
    fdesc.set_priors = exp_priors ;
    fdesc.NPARAMS    = 2 ;
    break ;
  case EXP_PLUSC : 
    fdesc.func       = fexp_plusc ;
    fdesc.F          = exp_plusc_f ;
    fdesc.dF         = exp_plusc_df ;
    fdesc.d2F        = exp_plusc_d2f ;
    fdesc.guesses    = exp_plusc_guesses ; 
    fdesc.set_priors = exp_plusc_priors ;
    fdesc.NPARAMS    = 3 ;
    break ;
  case PADE : 
    break ;
  case POLY : 
    break ;
  case SINH : 
    fdesc.func       = fsinh ;
    fdesc.F          = sinh_f ;
    fdesc.dF         = sinh_df ;
    fdesc.d2F        = sinh_d2f ;
    fdesc.guesses    = sinh_guesses ; 
    fdesc.set_priors = sinh_priors ;
    fdesc.NPARAMS    = 2 ;
    break ;
  }
  // allocate the fitfunction
  fdesc.f = allocate_ffunction( fdesc.NPARAMS , Ndata ) ;
  
  // corrfit
  fdesc.f.CORRFIT = CORRFIT ;

  return fdesc ;
}
