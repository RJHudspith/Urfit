/**
   @file resampled_ops.c
   @brief various mathematical operations on distributions
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fitfunc.h"

static void
NSAMPLES_ERR( const char *OP , const size_t a , const size_t b ) 
{
  printf( "[%s]\n" , OP ) ;
  printf( "the two distributions have different number of SAMPLES !!" ) ;
  printf( "%zu vs %zu \n" , a , b ) ;
  printf( "Leaving in disgust\n" ) ;
  return ;
}

static void
RESAMPLE_ERR( const char *OP , const size_t a , const size_t b ) 
{
  printf( "[%s]\n" , OP ) ;
  printf( "the two distributions have different sample types !!" ) ;
  printf( "%zu vs %zu \n" , a , b ) ;
  printf( "Leaving in disgust\n" ) ;
  return ;
}

// is always called
static void
resampled_checks( const char *OP , 
		  const struct resampled a , 
		  const struct resampled b ) {
  if( a.NSAMPLES != b.NSAMPLES ) { 
    NSAMPLES_ERR( OP , a.NSAMPLES , b.NSAMPLES ) ; 
  } else if( a.restype != b.restype ) {
    RESAMPLE_ERR( OP , a.restype , b.restype ) ; 
  }
  // and any other checks we might wish to perform ....
  return ;
}

// atomic, a += b for the distribution
void
add( struct resampled *a , 
     const struct resampled b ) 
{
  resampled_checks( "add" , *a , b ) ;
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] += b.resampled[i] ;
  }
  a -> avg += b.avg ;
  compute_err( a ) ;
  return ;
}

// atomic add by a (double)constant
void
add_constant( struct resampled *a , 
	      const double b ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] += b ;
  }
  a -> avg += b ;
  compute_err( a ) ;
  return ;
}

// atomic, a /= b for the distribution
void
divide( struct resampled *a , 
	const struct resampled b ) 
{
  resampled_checks( "divide" , *a , b ) ;
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] /= b.resampled[i] ;
  }
  a -> avg /= b.avg ;
  compute_err( a ) ;
  return ;
}

// atomic divide by a (double)constant
void
divide_constant( struct resampled *a , 
		 const double b ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] /= b ;
  }
  a -> avg /= b ;
  compute_err( a ) ;
  return ;
}

// equate two distributions
void
equate( struct resampled *a , 
	const struct resampled b )
{
  size_t l ;
  for( l = 0 ; l < b.NSAMPLES ; l++ ) {
    a -> resampled[l] = b.resampled[l] ;
  }
  a -> NSAMPLES = b.NSAMPLES ;
  a -> restype  = b.restype ;
  a -> avg = b.avg ;
  a -> err = b.err ;
  a -> err_hi = b.err_hi ;
  a -> err_lo = b.err_lo ;
  return ;
}

// equate a distribution to a constant
void
equate_constant( struct resampled *a ,
		 const double constant ,
		 const size_t NSAMPLES ,
		 const resample_type restype )
{
  size_t l ;
  for( l = 0 ; l < NSAMPLES ; l++ ) {
    a -> resampled[l] = constant ;
  }
  a -> NSAMPLES = NSAMPLES ;
  a -> restype  = restype ;
  a -> avg = constant ;
  a -> err = 0.0 ;
  a -> err_hi = constant ;
  a -> err_lo = constant ;
  return ;
}

// init to zero for (NULL,NSAMPLES,restype) else init to *d
struct resampled
init_dist( struct resampled *d , 
	   const int NSAMPLES , 
	   const resample_type restype )
{
  struct resampled sample ;
  sample.resampled = calloc( NSAMPLES , sizeof( double ) ) ;
  if( d == NULL ) {
    sample.NSAMPLES = NSAMPLES ;
    sample.restype  = restype ;
    sample.avg = 0.0 ;
    sample.err = 0.0 ;
    sample.err_hi = 0.0 ;
    sample.err_lo = 0.0 ;
  } else {
    equate( &sample , *d ) ;
  }
  return sample ;
}

// atomic, a *= b for the distribution
void
mult( struct resampled *a , 
      const struct resampled b ) 
{
  resampled_checks( "mult" , *a , b ) ;
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] *= b.resampled[i] ;
  }
  a -> avg *= b.avg ;
  compute_err( a ) ;
  return ;
}

// atomic multiply by a (double)constant
void
mult_constant( struct resampled *a , 
	       const double b ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] *= b ;
  }
  a -> avg *= b ;
  compute_err( a ) ;
  return ;
}

// atomic raise to a power a = a^b
void
raise( struct resampled *a , 
       const double b ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] = pow( a -> resampled[i] , b ) ;
  }
  a -> avg = pow( a -> avg , b ) ; 
  compute_err( a ) ;
  return ;
}

// atomic exponential a = res_acosh( a )
void
res_acosh( struct resampled *a ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] = acosh( a -> resampled[i] ) ;
  }
  a -> avg = acosh( a -> avg ) ;
  compute_err( a ) ;
  return ;
}

// atomic arcsinh a = res_asinh( a )
void
res_asinh( struct resampled *a ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] = asinh( a -> resampled[i] ) ;
  }
  a -> avg = asinh( a -> avg ) ;
  compute_err( a ) ;
  return ;
}

// atomic arctanh a = res_atanh( a )
void
res_atanh( struct resampled *a ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] = atanh( a -> resampled[i] ) ;
  }
  a -> avg = atanh( a -> avg ) ;
  compute_err( a ) ;
  return ;
}

// atomic exponential a = res_exp( a )
void
res_exp( struct resampled *a ) 
{
  int i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] = exp( a -> resampled[i] ) ;
  }
  a -> avg = exp( a -> avg ) ;
  compute_err( a ) ;
  return ;
}

// atomic logarithm a = log( a )
void
res_log( struct resampled *a ) 
{
  int i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] = log( a -> resampled[i] ) ;
  }
  a -> avg = log( a -> avg ) ;
  compute_err( a ) ;
  return ;
}

// atomic square root a = sqrt( a )
void
root( struct resampled *a ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] = sqrt( a -> resampled[i] ) ;
  }
  a -> avg = sqrt( a -> avg ) ;
  compute_err( a ) ;
  return ;
}

// atomic, a -= b for the distribution
void
subtract( struct resampled *a , 
	  const struct resampled b ) 
{
  resampled_checks( "subtract" , *a , b ) ;
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] -= b.resampled[i] ;
  }
  a -> avg -= b.avg ;
  compute_err( a ) ;
  return ;
}

// atomic subtract by a (double)constant
void
subtract_constant( struct resampled *a , 
		   const double b ) 
{
  size_t i ;
  for( i = 0 ; i < a -> NSAMPLES ; i++ ) {
    a -> resampled[i] -= b ;
  }
  a -> avg -= b ;
  compute_err( a ) ;
  return ;
}
