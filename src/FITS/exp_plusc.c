/**
   exponential Model f = ( A * exp(-lambda * i) + b - y_i )
 */
#include "gens.h"

double
fexp_plusc( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[1] * exp( -fparams[0] * X.X ) + fparams[2] ;
}

void
exp_plusc_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i ;
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA-> LT } ;
    f[i] = fexp_plusc( X , fparams , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
exp_plusc_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e = exp(-fparams[0] * t) ;
    df[0][i] = -t * fparams[1] * e ;
    df[1][i] = e ;
    df[2][i] = 1.0 ;
  }
  return ;
}

// second derivatives
void
exp_plusc_d2f( double **d2f , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  const double A = fparams[1] ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e = exp(-fparams[0] * t) ;
    d2f[0][i] = t*t*A*e ; d2f[1][i] = -t * e  ; d2f[2][i] = 0.0 ;
    d2f[3][i] = -t * e  ; d2f[4][i] = 0.0     ; d2f[5][i] = 0.0 ;
    d2f[6][i] = 0.0     ; d2f[7][i] = 0.0     ; d2f[8][i] = 0.0 ;
  }
  return ;
}

void
exp_plusc_guesses( double *fparams )
{
  if( fparams[0] == UNINIT_FLAG || fparams[1] == UNINIT_FLAG ||
      fparams[2] == UNINIT_FLAG ) {
  fparams[0] = 0.12 ; fparams[1] = 5.5 ; fparams[2] = 0.9 ;
  }
  return ;
}

void
exp_plusc_priors( double *priors , double *err_priors )
{
  /*
  priors[0] = 0.1 ; err_priors[1] = 0.005 ;
  priors[1] = 5.0 ; err_priors[0] = 0.1 ;
  priors[2] = 1.0 ; err_priors[2] = 0.05 ;
  */
}
