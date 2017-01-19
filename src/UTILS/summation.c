/**
   @file summation.c
   @brief round-off resistant summation routines
 */
#include "gens.h"

// perform a round-off resistant summation
double
kahan_summation( const double *data ,
		 const size_t Ndata )
{
  register double sum = 0.0 , t , c = 0.0 , tmp ;
  size_t i ;
  for( i = 0 ; i < Ndata ; i++ ) {
    tmp = data[i] - c ;
    t = sum + tmp ;
    c = ( t - sum ) - tmp ;
    sum = t ;
  }
  return sum ;
}

// computes the average and the variance in place
double
knuth_average( double *err ,
	       const double *data ,
	       const size_t Ndata )
{
  const double *y = data ;
  register double ave = 0.0 , delta ;
  *err = 0.0 ;
  size_t i ;
  for( i = 0 ; i < Ndata ; i++ ) {
    delta = *y - ave ;
    ave  += delta / ((double)i + 1.0 ) ;
    *err += delta * ( *y - ave ) ;
    y++ ;
  }
  return ave ;
}
