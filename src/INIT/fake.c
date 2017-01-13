/**
   @file fake_exp.c
   @brief generate some fake noisy data to test the fit codes
 */
#include "gens.h"

#include "ffunction.h"
#include "fit_chooser.h"
#include "pmap.h"
#include "stats.h"

#define xnoise (gsl_ran_gaussian( r , xsigma ) )
#define ynoise (gsl_ran_gaussian( r , ysigma ) )

// assumes x and y have been set
int
generate_fake_data( struct data_info *Data ,
		    const struct fit_info Fit ,
		    const double xsigma ,
		    const double ysigma )
{  
  // set up the gsl rng
  gsl_rng_env_setup( ) ;
  gsl_rng *r = gsl_rng_alloc( gsl_rng_default ) ;
  
  size_t i , j , k ;

  // allocate the pointers we are passing by reference
  Data -> x = malloc( Data -> Ntot * sizeof( struct resampled ) ) ;
  Data -> y = malloc( Data -> Ntot * sizeof( struct resampled ) ) ;
  
  // initialise the fit so we can get at the fit function
  struct fit_descriptor fdesc = init_fit( *Data , Fit ) ;

  // set the guesses
  fdesc.guesses( fdesc.f.fparams , fdesc.Nlogic ) ;
  
  size_t shift = 0 , p ;
  for( i = 0 ; i < Data -> Nsim ; i++ ) {
    for( j = shift ; j < shift + Data -> Ndata[i] ; j++ ) {

      Data -> x[j].resampled = malloc( Data -> Nboots * sizeof( double ) ) ;
      Data -> x[j].restype   = BootStrap ;
      Data -> x[j].NSAMPLES  = Data -> Nboots ;

      Data -> y[j].resampled = malloc( Data -> Nboots * sizeof( double ) ) ;
      Data -> y[j].restype   = BootStrap ;
      Data -> y[j].NSAMPLES  = Data -> Nboots ;

      // is a random "x" value
      const double x_prime = gsl_rng_uniform( r ) * Data -> Ndata[i] ;

      // fill the boots with noise
      for( k = 0 ; k < Data -> Nboots ; k++ ) {
	const double xx = x_prime + xnoise ;
        Data -> x[j].resampled[k] = xx ;
	// compute the fit func values
	double fparams[ fdesc.Nparam ] ;
	for( p = 0 ; p < fdesc.Nparam ; p++ ) {
	  fparams[ p ] = fdesc.f.fparams[ Fit.map[shift].p[p] ] ;
	}
	struct x_desc xdesc = { xx , Data -> LT } ;
	// evaluate the fit function and add some y-noise to it as well
	const double y_prime = fdesc.func( xdesc , fparams , fdesc.Nparam ) ;
        Data -> y[j].resampled[k] = y_prime * ( 1 + ynoise ) ;
      }

      compute_err( &(Data -> x[j]) ) ;
      compute_err( &(Data -> y[j]) ) ;
      printf( "%g %g || %g %g \n" , Data -> x[j].avg , Data -> x[j].err , Data -> y[j].avg , Data -> y[j].err ) ;
    }
    shift += Data -> Ndata[i] ;
  }

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  // free the rng
  gsl_rng_free( r ) ;
  
  return SUCCESS ;
}

#undef xnoise
#undef ynoise
