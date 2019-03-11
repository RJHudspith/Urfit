/**
   @file raw.c
   @brief raw data resampling
 */
#include "gens.h"
#include "summation.h"

// compute the average and returns the unnormalised variance
// using Knuth averaging method
void
average( double *ave , double *err , 
	 const double *meas , const size_t N )
{
  *ave = knuth_average( err , meas , N ) ;
  return ;
}

// compute the standard error and average
void
raw_err( struct resampled *replicas )
{
  average( &(replicas -> avg) ,
	   &(replicas -> err) ,
	   replicas -> resampled ,
	   replicas -> NSAMPLES ) ;
  replicas -> err = sqrt( replicas -> err / pow( replicas -> NSAMPLES - 1 , 2 ) ) ;
  replicas -> err_hi = replicas -> avg + replicas -> err ;
  replicas -> err_lo = replicas -> avg - replicas -> err ;
  return ;
}
