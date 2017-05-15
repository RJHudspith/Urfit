/**
   @file NR.c
   @brief Newton-Raphson method to solve f(x) = 0 for x from fit parameters
 */
#include "gens.h"

#include "fit_chooser.h"
#include "resampled_ops.h"
#include "stats.h"

// solve for zero using newton's method with finite difference derivative eval
// for parameters in one dimension
struct resampled
fit_zero( const struct resampled *Fit ,
	  const struct input_params *Input ,
	  const double initial_guess )
{
  const size_t Max_iters = 50 ;
  const double Precision = 1E-8 ;
  const double Hinit = 1E-1 , Hstep = 0.5 ;

  fprintf( stdout , "[NR] convergence criteria eps %e :: Max iters %zu \n" ,
	   Precision , Max_iters ) ;
  
  // init zero
  struct resampled zero = init_dist( NULL , Fit[0].NSAMPLES , Fit[0].restype ) ;
  
  // set up the fit again
  struct fit_descriptor fdesc = init_fit( Input -> Data , Input -> Fit ) ;

  size_t j , p , shift = 0 ;
  for( j = 0 ; j < Fit[0].NSAMPLES ; j++ ) {

    double NTOL = 1 , NTOL_PREV = UNINIT_FLAG , GUESS = initial_guess ;
    size_t Niter = 0 ;
    
    while( NTOL > Precision ) {
      
      struct x_desc xdesc = { GUESS , Input -> Data.LT[shift] ,
			      Input -> Fit.N , Input -> Fit.M } ;
      
      double fparams[ fdesc.Nparam ] ;

      // evaluate the fitfunc
      for( p = 0 ; p < fdesc.Nparam ; p++ ) {
	fparams[ p ] = Fit[p].resampled[j] ;
      }

      // numerically evaluate the derivative
      double DTOL = 1.0 , DPREV = UNINIT_FLAG , h = Hinit , der = 1.0 ;
      size_t niter = 0 ;
      
      while( DTOL > Precision ) {

	// five point stencil computation
	const struct x_desc xp2h = { GUESS+2*h , Input -> Data.LT[shift] ,
				     Input -> Fit.N , Input -> Fit.M } ;
	const struct x_desc xph  = { GUESS+h , Input -> Data.LT[shift] ,
				     Input -> Fit.N , Input -> Fit.M } ;
	const struct x_desc xmh  = { GUESS-h , Input -> Data.LT[shift] ,
				     Input -> Fit.N , Input -> Fit.M } ;
	const struct x_desc xm2h = { GUESS-2*h , Input -> Data.LT[shift] ,
				     Input -> Fit.N , Input -> Fit.M } ;

	const double trial = ( -fdesc.func( xp2h  , fparams , 0 ) + 
			       8*( fdesc.func( xph  , fparams , 0 ) -
				   fdesc.func( xmh , fparams , 0 ) )
			       +fdesc.func( xm2h  , fparams , 0 )
			       ) / ( 12*h ) ;

	DTOL = fabs( der - trial ) ;

	if( niter > Max_iters ) {
	  fprintf( stderr , "[NR] Inner der failing to converge %zu -> %e \n" ,
		   niter , DTOL ) ;
	  break ;
	}

	if( DPREV < DTOL ) {
	  fprintf( stderr , "[NR] Maximimum inner der convergence reached %e \n" ,
		   DTOL ) ;
	  break ;
	}

	h *= Hstep ;
	der = trial ;
	DPREV = DTOL ;
	niter++ ;
      }

      const double correction = 
	fdesc.func( xdesc , fparams , 0 ) / der ;
    
      NTOL = fabs( correction ) ;

      if( Niter > Max_iters ) {
	fprintf( stderr , "[NR] Newton-Raphson not converging after %zu iterations %e\n" ,
		 Niter , NTOL ) ;
	break ;
      }

      if( NTOL > NTOL_PREV ) {
	fprintf( stderr , "[NR] Maximum Newton-Raphson convergence reached %e\n" ,
		 NTOL ) ;
	break ;
      }
      
      GUESS -= correction ;

      #ifdef VERBOSE
      fprintf( stdout , "GUESS :: %f -> %e \n" , GUESS , NTOL ) ;
      #endif
      
      NTOL_PREV = NTOL ;
      Niter++ ;
    }
    zero.resampled[j] = GUESS ;
  }
  compute_err( &zero ) ;
  
  return zero ;
}
