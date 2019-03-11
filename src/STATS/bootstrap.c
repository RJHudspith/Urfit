/**
   @file bootstrap.c
   @brief bootstrapping code
 */
#include "gens.h"

#include "rng.h"        // for the bootstraps
#include "summation.h"

// perhaps do a naive bias correction -> I don't trust this!!
//#define BCORRECT

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

// qsort comparison function for the bootstrap
int 
comp_size_t( const void *elem1 , 
	     const void *elem2 ) 
{

  const size_t f = *( (size_t*)elem1 ) ;
  const size_t s = *( (size_t*)elem2 ) ;
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

  // compute the boot average and subtract the bias so that
  // it is the same as the true average
#ifdef BCORRECT
  double err ;
  const double bias = ( knuth_average( &err , replicas -> resampled ,
				       replicas -> NSAMPLES )
			- replicas -> avg ) ;
  size_t i ;
  for( i = 0 ; i < replicas -> NSAMPLES ; i++ ) {
    replicas -> resampled[i] -= bias ;
  }
#endif
  
  double *sorted = malloc( replicas -> NSAMPLES * sizeof( double ) ) ;
  memcpy( sorted , replicas -> resampled , sizeof( double ) * ( replicas -> NSAMPLES ) ) ;
  
  // sort the bootstraps
  qsort( sorted , replicas -> NSAMPLES , sizeof( double ) , comp ) ;

  // confidence bounds are at 1 sigma
  const double confidence = 68.2689492 ;
  const double omitted = 0.5 * ( 100. - confidence ) ;
  const size_t bottom = (size_t)( ( omitted * replicas -> NSAMPLES ) / 100. ) ;
  const size_t top = ( ( replicas -> NSAMPLES ) - 1 - bottom ) ;

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
  size_t i , shift = 0 ;

  // do a check to see if NSAMPLES is the same for all
  bool Nsampflag = true ;
  size_t bootmax = 0 , NSAMPLES = Input -> Data.x[0].NSAMPLES ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t j ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      if( Input -> Data.x[j].NSAMPLES > bootmax ) {
	bootmax = Input -> Data.x[j].NSAMPLES ;
      }
      if( Input -> Data.x[j].NSAMPLES != NSAMPLES ||
	  Input -> Data.y[j].NSAMPLES != NSAMPLES ) {
	Nsampflag = false ;
      }
    }
    shift = j ;
  }

  // perform a bootstrap
  double *xstrap = malloc( Input -> Data.Nboots * sizeof( double ) ) ;
  double *ystrap = malloc( Input -> Data.Nboots * sizeof( double ) ) ;

  // precompute rng sequence
  double *rng = NULL ;
  uint32_t **rng_idx = NULL ;
  if( Nsampflag == true ) {
    rng_idx = malloc( Input -> Data.Nboots * sizeof( uint32_t* ) ) ;
    for( i = 0 ; i < Input -> Data.Nboots ; i++ ) {
      rng_idx[i] = malloc( bootmax * sizeof( uint32_t ) ) ;
    }
  } else {
    rng = malloc( Input -> Data.Nboots * bootmax * sizeof( double ) ) ;
  }

  // precompute the random selections and sort them to help the cache
  rng_reseed() ;
  if( Nsampflag == true ) {
    size_t j ;
    for( i = 0 ; i < Input -> Data.Nboots ; i++ ) {
      for( j = 0 ; j < bootmax ; j++ ) {
	rng_idx[i][j] = (uint32_t)( rng_double() * bootmax ) ;
      }
      qsort( rng_idx[i] , bootmax , sizeof(uint32_t) , comp_size_t ) ;
    }
  } else {
    for( i = 0 ; i < Input -> Data.Nboots * bootmax ; i++ ) {
      rng[i] = rng_double() ;
    }
  }
  printf( "[BOOT] rng setup bootmax :: %zu \n" , bootmax ) ;
  
  // loop data doing the bootstrapping cannot be done in parallel yet
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    
    const double *x = Input -> Data.x[i].resampled ;
    const double *y = Input -> Data.y[i].resampled ;
    const size_t N = Input -> Data.x[i].NSAMPLES ;
    size_t k , l ;
    
    if( Nsampflag == true ) {
      
      for( k = 0 ; k < Input -> Data.Nboots ; k++ ) {
	
	// loop boots
	uint32_t *p = rng_idx[k] ;
	
	// loop raw data
	register double xsum = 0.0 , ysum = 0.0 ;
	for( l = 0 ; l < N ; l++ ) {
	  xsum += x[ *p ] ;
	  ysum += y[ *p ] ;
	  p++ ;
	}
	xstrap[k] = xsum / N ;
	ystrap[k] = ysum / N ;
      }
      
    } else {

      // loop boots
      double *p = rng ;
      
      for( k = 0 ; k < Input -> Data.Nboots ; k++ ) {
	
	// set the rng idx
	register size_t rng_idx = 0 ;
	
	// loop raw data
	register double xsum = 0.0 , ysum = 0.0 ;
	for( l = 0 ; l < N ; l++ ) {
	  rng_idx = (size_t)( ( *p ) * N ) ;
	  xsum += x[ rng_idx ] ;
	  ysum += y[ rng_idx ] ;
	  p++ ;
	}
	xstrap[k] = xsum / N ;
	ystrap[k] = ysum / N ;
      }
    }

    // reallocate and copy over
    Input -> Data.x[i].resampled =		\
      realloc( Input -> Data.x[i].resampled ,
	       Input -> Data.Nboots * sizeof( double ) ) ;
    Input -> Data.y[i].resampled =		\
      realloc( Input -> Data.y[i].resampled ,
	       Input -> Data.Nboots * sizeof( double ) ) ;

    for( k = 0 ; k < Input -> Data.Nboots ; k++ ) {
      Input -> Data.x[i].resampled[k] = xstrap[k] ;
      Input -> Data.y[i].resampled[k] = ystrap[k] ;
    }
    Input -> Data.x[i].restype = Input -> Data.y[i].restype = BootStrap ;
    Input -> Data.x[i].NSAMPLES = Input -> Data.Nboots ;
    Input -> Data.y[i].NSAMPLES = Input -> Data.Nboots ;
      
    // compute the error
    bootstrap_error( &(Input -> Data.x[i]) ) ;
    bootstrap_error( &(Input -> Data.y[i]) ) ;

    #ifdef VERBOSE
    printf( "BOOT %f %f || %f %f \n" ,
	    Input -> Data.x[i].avg , Input -> Data.x[i].err ,
	    Input -> Data.y[i].avg , Input -> Data.y[i].err ) ;
    #endif
  }
  
  // free the RNG sequence
  if( rng != NULL ) {
    free( rng ) ;
  }
  if( rng_idx != NULL ) {
    for( i = 0 ; i < Input -> Data.Nboots ; i++ ) {
      free( rng_idx[i] ) ;
    }
    free( rng_idx ) ;
  }
  
  free( xstrap ) ;
  free( ystrap ) ;

  return ;
}
