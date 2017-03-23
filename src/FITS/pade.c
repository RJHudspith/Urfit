/**
   PADE Model f = ( b + b*x + ... + b_n*x^n) / ( 1 + c*x + ... + c_mx^m )
 */
#include "gens.h"

#include "gls_bootfit.h"
#include "pade_coefficients.h"

double
fpade( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  double num = 0.0 , den = 0.0 ;
  for( i = 0 ; i < X.N ; i++ ) {
    num += fparams[i] * pow( X.X , i ) ;
  }
  for( i = 0 ; i < X.M ; i++ ) {
    den += fparams[X.N+i] * pow( X.X , i+1 ) ;
  }
  return num / ( 1.0 + den ) ;
}

void
pade_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  double par[ DATA -> Npars ] ;
  size_t i , p ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;
    for( p = 0 ; p < DATA -> Npars ; p++ ) {
      par[ p ] = fparams[ DATA -> map[i].p[p] ] ;
    }
    f[i] = fpade( X , par , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
pade_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double X = DATA -> x[i] ;
    double num = 0.0 , den = 0.0 ;
    for( j = 0 ; j < DATA -> N ; j++ ) {
      num += fparams[ DATA -> map[i].p[j] ] * pow( X , j ) ;
    }
    for( j = 0 ; j < DATA -> M ; j++ ) {
      den += fparams[ DATA -> map[i].p[ DATA -> N +j] ] * pow( X , j+1 ) ;
    }
    den += 1.0 ;
    const double fac = num / den ;
    // numerator derivatives
    for( j = 0 ; j < DATA -> N ; j++ ) {
      df[j][i] = pow( X , j ) / den ;
    }
    for( j = 0 ; j < DATA -> M ; j++ ) {
      df[j+DATA -> N][i] = -pow( X , j+1 ) * fac / den ;
    }
  }
  return ;
}

// second derivatives
void
pade_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
pade_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit ) 
{
  fparams[0] = 1.0 ;
  fparams[1] = 0.0 ;
  fparams[2] = -0.5 ;
  /*
  fparams[0] = -0.18 ;
  fparams[1] = -0.08 ;
  fparams[2] = 0.0007 ;
  fparams[3] = 0.80 ;
  fparams[4] = 0.019 ;
  fparams[5] = -0.00014 ;
  */

  return ;
}
