/**
   @file bootstrap.c
   @brief bootstrapping code
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fitfunc.h"
#include "rng.h"      // for the bootstraps
#include "stats.h"

// either we boot average or use the ensemble average ...
#define BOOT_AVERAGE

static double *boot_order ;

// qsort comparison function for the bootstrap
static int 
comp( const void *elem1 , 
      const void *elem2 ) 
{
  const double f = *( (double*)elem1 ) ;
  const double s = *( (double*)elem2 ) ;
  if (f > s) return  1 ;
  if (f < s) return -1 ;
  return 0 ;
}

// free the bootstrap order
void
free_boot_order( void )
{
  free( boot_order ) ;
  return ;
}

// initialise the bootstrap order
void
init_boot_order( const size_t NBOOTS , 
		 const size_t NRAW )
{
  rng_reseed( ) ;
  boot_order = (double*)malloc( NBOOTS * NRAW * sizeof( double ) ) ;
  size_t i ;
  for( i = 0 ; i < NBOOTS * NRAW  ; i++ ) {
    boot_order[ i ] = rng_double( ) ;
  }
  return ;
}

// compute the error for the bootstrap assumes the data has
// been bootstrapped ...
void
bootstrap_error( struct resampled *replicas )
{
  double sorted[ replicas -> NSAMPLES ] ;
  memcpy( sorted , replicas -> resampled , sizeof( double ) * ( replicas -> NSAMPLES ) ) ;

  // sort the bootstraps
  qsort( sorted , replicas -> NSAMPLES , sizeof( double ) , comp ) ;

  // average is better behaved when sorted ...
  double ave = 0.0 ;
  size_t i ;
  for( i = 0 ; i < ( replicas -> NSAMPLES ) ; i++ ) {
    ave += sorted[ i ] ;
  }
  ave /= (double)( replicas -> NSAMPLES ) ;

  // confidence bounds are at 1 sigma
  const double confidence = 68.2689492 ;
  const double omitted = 0.5 * ( 100. - confidence ) ;
  const size_t bottom = (size_t)( ( omitted * replicas -> NSAMPLES ) / 100. ) ;
  const size_t top = ( ( replicas -> NSAMPLES ) - 1 - bottom ) ;

#ifdef BOOT_AVERAGE
  replicas -> avg = ave ;
#endif
  replicas -> err_hi = sorted[ top ] ;
  replicas -> err_lo = sorted[ bottom ] ;
  // symmetrized error
  replicas -> err    = 0.5 * ( replicas -> err_hi - replicas -> err_lo ) ;
  return ;
}

// and the bootstrap analysis , assumes we have NBOOTS of space in Bootstrap
void
bootstrap( struct resampled *Bootstrap ,
	   const struct resampled Raw )
{
  // reseed the rng ...
  rng_reseed( ) ;
  // loop the bootstrap index
  size_t i , j ;
  for( i = 0 ; i < Bootstrap -> NSAMPLES ; i++ ) {
    Bootstrap -> resampled[ i ] = 0.0 ;    
    // loop the data
    for( j = 0 ; j < Raw.NSAMPLES ; j++ ) {
      Bootstrap -> resampled[ i ] += 
	Raw.resampled[ (size_t)( boot_order[ j + Raw.NSAMPLES * i ] * 
			      ( Raw.NSAMPLES - 1 ) ) ] ;
    }
    Bootstrap -> resampled[ i ] /= (double)Raw.NSAMPLES ;
  }
#ifndef BOOT_AVERAGE
  Bootstrap -> avg = Raw.avg ;
#endif
  bootstrap_error( Bootstrap ) ;
  return ;
}
