/**
   @file fake_exp.c
   @brief generate some fake noisy data to test the fit codes
 */
#include "gens.h"

#include "ffunction.h"
#include "fit_chooser.h"
#include "init.h"
#include "pmap.h"
#include "stats.h"
#include "resampled_ops.h"

struct resampled *
generate_fake_boot( const size_t N,
		    const size_t Nsamples,
		    const double yarr[N] ,
		    const double dyarr[N] )
{
  //printf( "%zu %zu\n" , N , Nsamples ) ;

  
  struct resampled *y = malloc( N * sizeof( struct resampled ) ) ;
  
  gsl_rng *r = NULL ;
  // set up the gsl rng
  gsl_rng_env_setup( ) ;
  r = gsl_rng_alloc( gsl_rng_default ) ;

  printf( "In fake boot\n" ) ;
  for( size_t i = 0 ; i < N ; i++ ) {
    y[i] = init_dist( NULL , Nsamples , BootStrap ) ;

    printf( "y[i] init\n" ) ;

    y[i].avg = yarr[i] ;
    for( size_t j = 0 ; j < Nsamples ; j++ ) {
      y[i].resampled[j] = yarr[i] + gsl_ran_gaussian( r , dyarr[i] ) ;
    }

    printf( "y[i] error\n" ) ;
    
    compute_err( &y[i] ) ;

    printf( "%zu done\n" , N ) ;
  }
  
  return y ;
}

// assumes x and y have been set
int
generate_fake_data( struct data_info *Data ,
		    struct fit_info Fit ,
		    struct traj *Traj ,
		    const double xsigma ,
		    const double ysigma )
{
  const size_t Nmeas = 200 ;
  const size_t Ndata = 30 ;
  
  size_t i , j , k ;
  Data -> Ndata = malloc( Data -> Nsim * sizeof( size_t ) ) ;
  Data -> Ntot = 0 ;
  for( i = 0 ; i < Data -> Nsim ; i++ ) {
    Data -> Ndata[i] = Ndata ;
    Data -> Ntot += Data -> Ndata[i] ;
  }

  gsl_rng *r = NULL ;

   // set Lt
  if( init_LT( Data , Traj ) == FAILURE ) {
    goto memfree ;
  }
  
  // set up the gsl rng
  gsl_rng_env_setup( ) ;
  r = gsl_rng_alloc( gsl_rng_default ) ;
  
  // allocate the pointers we are passing by reference
  Data -> x = malloc( Data -> Ntot * sizeof( struct resampled ) ) ;
  Data -> y = malloc( Data -> Ntot * sizeof( struct resampled ) ) ;
  
  // initialise the fit so we can get at the fit function
  struct fit_descriptor fdesc = init_fit( *Data , Fit ) ;

  // set the param map
  Fit.map = parammap( *Data , Fit ) ;

  // set the fit parameters
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    fdesc.f.fparams[i] = ( i+1 ) ;
    printf( "FAKED_%zu %f \n" , i , fdesc.f.fparams[i] ) ;
  }
  
  size_t shift = 0 , p ;
  for( i = 0 ; i < Data -> Nsim ; i++ ) {
    for( j = shift ; j < shift + Data -> Ndata[i] ; j++ ) {

      Data -> x[j].resampled = malloc( Nmeas * sizeof( double ) ) ;
      Data -> x[j].restype   = Raw ;
      Data -> x[j].NSAMPLES  = Nmeas ;

      Data -> y[j].resampled = malloc( Nmeas * sizeof( double ) ) ;
      Data -> y[j].restype   = Raw ;
      Data -> y[j].NSAMPLES  = Nmeas ;

      // is a random "x" value
      const double x_prime = j ; //gsl_rng_uniform( r ) * Ndata / 20 ;

      // fill the boots with noise
      for( k = 0 ; k < Nmeas ; k++ ) {
	const double xx = x_prime + gsl_ran_gaussian( r , xsigma ) ;
        Data -> x[j].resampled[k] = xx ;
	// compute the fit func values
	double fparams[ fdesc.Nparam ] ;
	for( p = 0 ; p < fdesc.Nparam ; p++ ) {
	  fparams[ p ] = fdesc.f.fparams[ Fit.map[shift].p[p] ] ;
	}
	struct x_desc xdesc = { xx , Data -> LT[j] , Fit.N , Fit.M } ;
	// evaluate the fit function and add some y-noise to it as well
	const double y_prime = fdesc.func( xdesc , fparams , fdesc.Nparam ) ;
        Data -> y[j].resampled[k] = y_prime *
	  ( 1 + gsl_ran_gaussian( r , ysigma ) ) ;
      }

      compute_err( &(Data -> x[j]) ) ;
      compute_err( &(Data -> y[j]) ) ;
      printf( "%g %g || %g %g \n" , Data -> x[j].avg , Data -> x[j].err , Data -> y[j].avg , Data -> y[j].err ) ;
    }
    shift += Data -> Ndata[i] ;
  }

 memfree :
  
  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  // free the rng
  gsl_rng_free( r ) ;
  
  // free the parameter map
  if( Fit.map != NULL ) {
    free_pmap( Fit.map , Data -> Ntot ) ;
  }
  
  return SUCCESS ;
}

#undef xnoise
#undef ynoise
