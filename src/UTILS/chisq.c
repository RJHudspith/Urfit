/**
   @file chisq.c
   @brief compute the chisq
 */
#include "gens.h"
#include "summation.h"

// compute the chisq using Kahan summation
double
compute_chisq( const struct ffunction f ,
	       const double **W ,
	       const corrtype CORRFIT )
{
  // pointers and stuff
  const double *pf = f.f , *pW = *W ;
  register double chisq = 0.0 ;
  size_t i , j , Nsum = 1 ;

  // allocate the array we want to sum
  switch( CORRFIT ) {
  case UNWEIGHTED : Nsum = f.N ; break ;
  case UNCORRELATED : Nsum = f.N ; break ;
  case CORRELATED : Nsum = f.N * f.N ; break ;
  }
  double y[ Nsum ] ;

  // compute the chisq within the switch
  switch( CORRFIT ) {
  case UNWEIGHTED :
    for( i = 0 ; i < f.N ; i++ ) {
      y[i] = *pf * ( *pf ) ;
      pf++ ;
    }
    break ;
  case UNCORRELATED :
    for( i = 0 ; i < f.N ; i++ ) {
      y[i] = *pf * (*pW) * ( *pf ) ;
      pf++ , pW++ ;
    }
    break ;
  case CORRELATED :
    for( i = 0 ; i < f.N ; i++ ) {
      register const double fi = f.f[i] ;
      // point to the data and the W matrix
      pf = f.f ; pW = *( W + i ) ;
      for( j = 0 ; j < f.N ; j++ ) {
	y[j+f.N*i] = fi * ( *pW ) * ( *pf ) ;
	pf++ , pW++ ;
      }
    }
    break ;
  }

  // perform the round off resistant summation
  chisq = kahan_summation( y , Nsum ) ;
  
  // add priors to the chisq
  for( i = 0 ; i < f.NPARAMS ; i++ ) {
    if( f.Prior[i].Initialised == true ) {
      const double fac = ( f.fparams[i] - f.Prior[i].Val ) / f.Prior[i].Err ;
      chisq += fac * fac ;
    }
  }
  return chisq ;
}
