/**
   @file pvalue.c
   @brief compute the p value from a statistical distribution
 */
#include "gens.h"

// return a p-value for our data given some "true value"
// S0, note that when S0 is the average of our data we
// necessarily get a p-value of 0
double
pvalue( struct resampled result ,
	const double S0 )
{
  const size_t N = result.NSAMPLES ;
  
  // create a copy
  double sorted[ N ] ;

  // recentre the array on the statistic S0
  for( i = 0 ; i < N ; i++ ) {
    sorted[ i ] = result.resampled[i] - result.avg + S0 ;
  }
  
  // sort the data smallest to largest
  qsort( sorted , N , sizeof( double ) , comp ) ;

  // compute the difference between the true mean and our observed mean
  const double diff = fabs( result.avg - S0 ) ;
  
  double var = 0.0 ;
  size_t i , p = 0 ;
  switch( result.restype ) {
  case BOOTSTRAP :
  case JACKKNIFE :
    // count how many values fall outside the difference
    for( i = 0 ; i < N ; i++ ) {
      if( ( sorted[ i ] < ( s0 - diff ) ) ||
	  ( sorted[ i ] > ( s0 + diff ) ) ) {
	p++ ;
      }
    }
    return ( 1 + p / ( (double)N + 1 ) ) ;
  default : return -1 ;
  }
  return -1 ;
}
