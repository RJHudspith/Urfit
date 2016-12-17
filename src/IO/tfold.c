/**
   @file tfold.c
   @brief time correlation function folding
 */
#include <complex.h>

#include "gens.h"

int
time_fold( struct resampled *sample ,
	   const double complex *C ,
	   const size_t LT ,
	   const foldtype fold ,
	   const size_t meas )
{
  const size_t L2 = LT/2 ;
  size_t t ;
  switch( fold ) {
  case PLUS_PLUS :
    sample[0].resampled[meas] = fabs( creal( C[0] ) ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * creal( C[t] + C[LT-t] ) ;
    }
    break ;
  case PLUS_MINUS :
    sample[0].resampled[meas] = fabs( creal( C[0] ) ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * creal( C[t] - C[LT-t] ) ;
    }
    break ;
  case MINUS_PLUS :
    sample[0].resampled[meas] = fabs( creal( C[0] ) ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * creal( -C[t] + C[LT-t] ) ;
    }
    break ;
  case MINUS_MINUS :
    sample[0].resampled[meas] = fabs( creal( C[0] ) ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * creal( -C[t] + C[LT-t] ) ;
    }
    break ;
  case NOFOLD :
    for( t = 0 ; t < LT ; t++ ) {
      sample[t].resampled[meas] = creal( C[t] ) ;
    }
    break ;
  }
  return SUCCESS ;
}
