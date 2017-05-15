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
#include "rng.h"
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
  replicas -> err = sqrt( replicas -> err / ( replicas -> NSAMPLES - 1 ) ) ;
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
int
resample_data( struct input_params *Input )
{
  size_t i ;
  bool must_resample = false ;
  
  // compute the error
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    compute_err( &Input -> Data.x[i] ) ;
    compute_err( &Input -> Data.y[i] ) ;
    if( Input -> Data.x[i].restype != Input -> Data.Restype ||
	Input -> Data.x[i].restype != Input -> Data.Restype )
      {
	must_resample = true ;
      }
  }

  // bootstrap it
  if( must_resample == true ) {
    switch( Input -> Data.Restype ) {
    case BootStrap :
      #ifdef VERBOSE
      fprintf( stdout , "[STATS] Bootstrapping \n" ) ;
      #endif
      init_rng( 123456 ) ;
      bootstrap_full( Input ) ;
      free_rng() ;
      break ;
    case JackKnife :
      #ifdef VERBOSE
      fprintf( stdout , "[STATS] JackKnifing \n" ) ;
      #endif
      jackknife_full( Input ) ;
      break ;
    case Raw :
      #ifdef VERBOSE
      fprintf( stdout , "[STATS] Raw data\n" ) ;
      #endif
      break ;
    }
  }

  #ifdef VERBOSE
  fprintf( stdout , "[STATS] Resampling finished\n" ) ;
  #endif

  return SUCCESS ;
}
