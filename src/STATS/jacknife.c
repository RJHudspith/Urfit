/**
   @file jacknife.c
   @brief jacknife resampling
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "gens.h"
#include "stats.h"

// the jackknife error definition incorporates the code for the average and variance
void
jackknife_error( struct resampled *replicas )
{
  // average and variance come from the usual defs
  double ave , err ;
  average( &ave , &err , replicas -> resampled , replicas -> NSAMPLES ) ;

  err = sqrt( err * ( 1.0 - 1.0 / (double)replicas -> NSAMPLES ) ) ;
  replicas -> avg = ave ;
  replicas -> err_hi = ave + err ;
  replicas -> err_lo = ave - err ;
  replicas -> err    = err ;
  return ;
}

// perform a jackknife on the Raw data
// assumes Jackknife has had NRAW - 1 allocation
void
jackknife( struct resampled *Jackknife ,
	   const struct resampled Raw )
{
  const size_t N = Raw.NSAMPLES ;
  Jackknife -> NSAMPLES = N ; // should be set anyway ..

  register double sum = 0.0 ;
  size_t i ;
  for( i = 0 ; i < N ; i++ ) {
    sum += Raw.resampled[i] ; // compute the sum
  }
  // subtract each element from the sum
  const double NORM = 1.0 / ( (double)N - 1.0 ) ;
  for( i = 0 ; i < N ; i++ ) {
    Jackknife -> resampled[ i ] = ( sum - Raw.resampled[i] ) * NORM ;
  }
  jackknife_error( Jackknife ) ;
  return ;
}
