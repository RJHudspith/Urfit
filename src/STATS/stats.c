/**
   @file stats.c
   @brief statistical sampling
 */
#include <stdlib.h>
#include <stdio.h>

#include "gens.h"
#include "bootstrap.h"
#include "jacknife.h"
#include "resampled_ops.h"
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

void
raw_err( struct resampled *replicas )
{
  average( &(replicas -> avg) ,
	   &(replicas -> err) ,
	   replicas -> resampled ,
	   replicas -> NSAMPLES ) ;
  replicas -> err = sqrt( replicas -> err ) / replicas -> NSAMPLES ;
  replicas -> err_hi = replicas -> avg + replicas -> err ;
  replicas -> err_lo = replicas -> avg - replicas -> err ;
  return ;
}

// bins raw data
struct resampled
bin_the_data( const struct resampled RAW ,
	      const size_t binning )
{
  const size_t NBINNED = (size_t)RAW.NSAMPLES / (size_t)binning ;
  struct resampled BINNED ;

  BINNED.resampled = malloc( NBINNED * sizeof( double ) ) ;
  BINNED.NSAMPLES = NBINNED ;
  BINNED.restype  = Raw;

  // OK, so we simply average within each bin
  size_t j , k ;
  for( j = 0 ; j < NBINNED ; j++ ) {
    // loop the elemnts within the bin, combining for an average
    register double bin_ave = 0.0 ;
    for( k = 0 ; k < binning ; k++ ) {
      bin_ave += RAW.resampled[ k + j*binning ] ;
    }
    BINNED.resampled[j] = bin_ave / (double)binning ;
  }
  return BINNED ;
}

// a wrapper for the resampling
void
compute_err( struct resampled *replicas )
{
  switch( replicas -> restype ) {
  case Raw :
    raw_err( replicas ) ;
    break ;
  case JackKnife :
    jackknife_error( replicas ) ;
    break ;
  case BootStrap :
    bootstrap_error( replicas ) ;
    break ;
  }
  return ;
}

// bootstrap or jackknife or whatever
struct resampled *
resample_data( const struct resampled *RAW ,
	       const size_t N ,
	       const resample_type restype ,
	       const int NBOOTS )
{
  return NULL ;
}
