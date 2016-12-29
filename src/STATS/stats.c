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

// compute the average and returns the unnormalised variance
// using Knuth averaging method
void
average( double *ave , double *err , 
	 const double *meas , const size_t N )
{
  size_t i ;
  *err = *ave = 0.0 ;
  for( i = 0 ; i < N ; i++ ) {
    const double delta = meas[i] - *ave ;
    *ave += delta / ((double)i+1.0) ;
    *err += delta * ( meas[i] - *ave ) ; 
  }
  return ;
}

// bins the raw data
struct resampled
bin_the_data( const struct resampled RAW ,
	      const size_t binning )
{
  const size_t NBINNED = (size_t)RAW.NSAMPLES / (size_t)binning ;
  struct resampled BINNED ;

  BINNED.resampled = malloc( NBINNED * sizeof( double ) ) ;
  BINNED.NSAMPLES = NBINNED ;
  BINNED.restype  = RAWDATA ;

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
  case RAWDATA :
    //raw_err( replicas ) ;
    break ;
  case JACKDATA :
    jackknife_error( replicas ) ;
    break ;
  case BOOTDATA :
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
  struct resampled *BOOT = malloc( N * sizeof( struct resampled ) ) ;
  const int NRAW = RAW[0].NSAMPLES ;
  size_t i ;

  printf( "[STATS] NBOOTS :: %d \n" , NBOOTS ) ;

  // and resample could be done in parallel because we have the rng
  // order as a look up table
  #pragma omp parallel for private(i)
  for( i = 0 ; i < N ; i++ ) {
    // resampling
    size_t k ;

    if( RAW[i].restype != RAWDATA ) {

      BOOT[i].resampled = (double*)malloc( RAW[i].NSAMPLES * sizeof( double ) ) ;
      equate( &BOOT[i] , RAW[i] ) ;

    } else {

      BOOT[i].restype = restype ;

      switch( restype ) {
      case RAWDATA :
	BOOT[i].resampled = (double*)malloc( NRAW * sizeof( double ) ) ;
	BOOT[i].NSAMPLES  = NRAW ;
	for( k = 0 ; k < NRAW ; k++ ) {
	  BOOT[ i ].resampled[ k ] = RAW[ i ].resampled[ k ] ;
	}
	//raw_err( &BOOT[i] ) ;
	break ;
      case BOOTDATA :
	BOOT[i].resampled = (double*)malloc( NBOOTS * sizeof( double ) ) ;
	BOOT[i].NSAMPLES  = NBOOTS ;
	bootstrap( &BOOT[i] , RAW[i] ) ;
	break ;
      case JACKDATA :
	// I manipulate NBOOTS to be NRAW-1 in the input file
	BOOT[i].resampled = (double*)malloc( NBOOTS * sizeof( double ) ) ;
	BOOT[i].NSAMPLES  = NBOOTS ;
	jackknife( &BOOT[i] , RAW[i] ) ;
	break ;
      }
      // break for the check on if it is rawdata or not
    }
  }
  return BOOT ;
}

// wraps the above into an array
struct resampled **
resample_array( const struct resampled **RAW ,
		const size_t NSLICES ,
		const size_t *NDATA ,
		const size_t resample , 
		const size_t NBOOTS )
{
  struct resampled **BOOT = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  size_t i ;
  for( i = 0 ; i < NSLICES ; i++ ) {
    BOOT[i] = resample_data( RAW[i] , 
			     NDATA[i] , 
			     resample , 
			     NBOOTS ) ;
  }
  return BOOT ;
}
