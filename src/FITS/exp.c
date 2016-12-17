/**
   exponential Model f = ( A * exp(-lambda * i) - y_i )
 */
#include "gens.h"

double
fexp( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[1] * exp( -fparams[0] * X.X ) ;
}

void
exp_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT } ;
    f[i] = fexp( X , fparams , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
exp_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e = exp(-fparams[0] * t ) ;
    df[0][i] = -t * fparams[1] * e ;
    df[1][i] = e ;
  }
  return ;
}

// second derivatives
void
exp_d2f( double **d2f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  const double A = fparams[1] ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e = exp(-fparams[0] * t) ;
    d2f[0][i] = t*t*A*e  ; d2f[1][i] = -t*e ; 
    d2f[2][i] = -t*e     ; d2f[3][i] = 0.0  ; 
  }
  return ;
}

void
exp_guesses( double *fparams )
{
  if( fparams[0] == UNINIT_FLAG && fparams[1] == UNINIT_FLAG ) {
    fparams[0] = 0.2 ; fparams[1] = 6.0 ;
  }
  return ;
}

void
exp_priors( double *priors , double *err_priors )
{
  return ;
}
