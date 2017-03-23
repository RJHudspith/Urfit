/**
   @file bootstrap.c
   @brief bootstrapping code
 */
#include "gens.h"

#include "rng.h"        // for the bootstraps
#include "stats.h"
#include "summation.h"

// either we boot average or use the ensemble average ...
//#define BOOT_AVERAGE

// qsort comparison function for the bootstrap
int 
comp( const void *elem1 , 
      const void *elem2 ) 
{

  const double f = *( (double*)elem1 ) ;
  const double s = *( (double*)elem2 ) ;
  if (f > s) { return  1 ; }
  if (f < s) { return -1 ; }
  return 0 ;
}

// compute the error for the bootstrap assumes the data has
// been bootstrapped ...
void
bootstrap_error( struct resampled *replicas )
{
  if( replicas -> NSAMPLES == 0 ) return ;
  
  double *sorted = malloc( replicas -> NSAMPLES * sizeof( double ) ) ;
  memcpy( sorted , replicas -> resampled , sizeof( double ) * ( replicas -> NSAMPLES ) ) ;
  
  // sort the bootstraps
  qsort( sorted , replicas -> NSAMPLES , sizeof( double ) , comp ) ;

  // confidence bounds are at 1 sigma
  const double confidence = 68.2689492 ;
  const double omitted = 0.5 * ( 100. - confidence ) ;
  const size_t bottom = (size_t)( ( omitted * replicas -> NSAMPLES ) / 100. ) ;
  const size_t top = ( ( replicas -> NSAMPLES ) - 1 - bottom ) ;

#ifdef BOOT_AVERAGE
  // is the middle of the bootstrap distribution ok?
  replicas -> avg = sorted[ replicas -> NSAMPLES / 2 ] ;
#endif
  replicas -> err_hi = sorted[ top ] ;
  replicas -> err_lo = sorted[ bottom ] ;
  // symmetrized error
  replicas -> err    = 0.5 * ( replicas -> err_hi - replicas -> err_lo ) ;

  free( sorted ) ;
  
  return ;
}

// perform a bootstrap resampling on the data
void
bootstrap_full( struct input_params *Input )
{
  size_t i , j = 0 , k , l , shift = 0 ;

  // perform a bootstrap
  double *xstrap = malloc( Input -> Data.Nboots * sizeof( double ) ) ;
  double *ystrap = malloc( Input -> Data.Nboots * sizeof( double ) ) ;

  // get the maximum
  size_t bootmax = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      if( Input -> Data.x[j].NSAMPLES > bootmax ) {
	bootmax = Input -> Data.x[j].NSAMPLES ;
      }
    }
    shift = j ;
  }

  // precompute rng sequence
  double *rng = malloc( Input -> Data.Nboots * bootmax * sizeof( double* ) ) ;
  rng_reseed() ;
  for( i = 0 ; i < Input -> Data.Nboots * bootmax ; i++ ) {
    rng[i] = rng_double() ;
  }
  printf( "[BOOT] rng setup bootmax :: %zu \n" , bootmax ) ;

  // loop data doing the bootstrapping
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      const size_t N = Input -> Data.x[j].NSAMPLES ;	
      double *x = malloc( N * sizeof( double ) ) ;
      double *y = malloc( N * sizeof( double ) ) ;

      // loop boots
      double *p = rng ;
      
      for( k = 0 ; k < Input -> Data.Nboots ; k++ ) {

	// set the rng idx
	register size_t rng_idx = 0 ;

	// loop raw data
	for( l = 0 ; l < N ; l++ ) {
	  rng_idx = (size_t)( ( *p ) * N ) ;
	  x[l] = Input -> Data.x[j].resampled[ rng_idx ] ;
	  y[l] = Input -> Data.y[j].resampled[ rng_idx ] ;
	  p++ ;
	}
	
	xstrap[k] = kahan_summation( x , N ) / N ;
	ystrap[k] = kahan_summation( y , N ) / N ;	
      }

      // reallocate and copy over
      Input -> Data.x[j].resampled = \
	realloc( Input -> Data.x[j].resampled ,
		 Input -> Data.Nboots * sizeof( double ) ) ;
      Input -> Data.y[j].resampled = \
	realloc( Input -> Data.y[j].resampled ,
		 Input -> Data.Nboots * sizeof( double ) ) ;

      for( k = 0 ; k < Input -> Data.Nboots ; k++ ) {
	Input -> Data.x[j].resampled[k] = xstrap[k] ;
	Input -> Data.y[j].resampled[k] = ystrap[k] ;
      }
      Input -> Data.x[j].restype = Input -> Data.y[j].restype = BootStrap ;
      Input -> Data.x[j].NSAMPLES = Input -> Data.Nboots ;
      Input -> Data.y[j].NSAMPLES = Input -> Data.Nboots ;

      #ifndef BOOT_AVERAGE
      x = realloc( x , Input -> Data.Nboots * sizeof( double ) ) ;
      y = realloc( y , Input -> Data.Nboots * sizeof( double ) ) ;
      for( l = 0 ; l < Input -> Data.Nboots ; l++ ) {
	x[l] = Input -> Data.x[j].resampled[l] ;
	y[l] = Input -> Data.y[j].resampled[l] ;
      }
      Input -> Data.x[j].avg = kahan_summation( x , Input -> Data.Nboots ) / Input -> Data.Nboots ;
      Input -> Data.y[j].avg = kahan_summation( y , Input -> Data.Nboots ) / Input -> Data.Nboots ;
      #endif

      // free the temp vectors
      free( x ) ; free( y ) ;
      
      // compute the error
      bootstrap_error( &(Input -> Data.x[j]) ) ;
      bootstrap_error( &(Input -> Data.y[j]) ) ;

      #ifdef VERBOSE
      printf( "BOOT %f %f || %f %f \n" ,
	      Input -> Data.x[j].avg , Input -> Data.x[j].err ,
	      Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
      #endif
    }
    shift += Input -> Data.Ndata[i] ;
  }
  
  // free the RNG sequence
  free( rng ) ;
  
  free( xstrap ) ;
  free( ystrap ) ;

  return ;
}
