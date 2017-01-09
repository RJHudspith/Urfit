#include "gens.h"

// macro expansion for the kahan summation
#define KAHAN(chisq,y,t,c) {				\
    t = chisq + y ; c = ( t - chisq ) - y ; chisq = t ; \
  }

// compute the chisq using Kahan summation
double
compute_chisq( struct ffunction f ,
	       const double **W ,
	       const corrtype CORRFIT )
{
  // Kahan summation parameters
  register double chisq = 0.0 , c = 0.0 , y , t ;
  size_t i , j ;

  const double *pf = f.f , *pW = *W ;
  
  switch( CORRFIT ) {
  case UNWEIGHTED :
    for( i = 0 ; i < f.N ; i++ ) {
      y = *pf * ( *pf ) - c ;
      KAHAN(chisq,y,t,c) ;
      pf++ ;
    }
    break ;
  case UNCORRELATED :
    for( i = 0 ; i < f.N ; i++ ) {
      y = *pf * (*pW) * ( *pf ) - c ;
      KAHAN(chisq,y,t,c) ;
      pf++ , pW++ ;
    }
    break ;
  case CORRELATED :
    for( i = 0 ; i < f.N ; i++ ) {
      const register double fi = f.f[i] ;
      // point to the data and the W matrix
      pf = f.f ; pW = *( W + i ) ;
      for( j = 0 ; j < f.N ; j++ ) {
	y = fi * ( *pW ) * ( *pf ) - c ;
	KAHAN(chisq,y,t,c) ;
	pf++ , pW++ ;
      }
    }
    break ;
  }
  
  // add priors to the chisq
  for( i = 0 ; i < f.NPARAMS ; i++ ) {
    if( f.prior[i] != UNINIT_FLAG ) {
      const double fac = ( f.fparams[i] - f.prior[i] ) / f.err_prior[i] ;
      chisq += fac * fac ;
    }
  }
  return (double)chisq ;
}

#undef KAHAN
