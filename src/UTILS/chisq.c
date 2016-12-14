#include "gens.h"

double
compute_chisq( struct ffunction f ,
	       const double **W ,
	       const corrtype CORRFIT )
{
  register double chisq = 0.0 ;
  size_t i , j ;
  for( i = 0 ; i < f.N ; i++ ) {
    switch( CORRFIT ) {
    case UNWEIGHTED : 
     chisq += f.f[i] * f.f[i] ;
      break ;
    case UNCORRELATED :
      chisq += f.f[i] * W[i][i] * f.f[i] ;
      break ;
    case CORRELATED :
      for( j = 0 ; j < f.N ; j++ ) {
	chisq += f.f[i] * W[i][j] * f.f[j] ;
      }
      break ;
    }
  }
  // add priors to the chisq
  for( i = 0 ; i < f.NPARAMS ; i++ ) {
    if( f.prior[i] != UNINIT_FLAG ) {
      const double fac = ( f.fparams[i] - f.prior[i] ) / f.err_prior[i] ;
      chisq += fac * fac ;
    }
  }
  return chisq ;
}
