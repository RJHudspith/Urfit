/**
   exponential Model f = ( A * exp(-lambda * i) + b - y_i )
 */
#include "gens.h"

void
exp_f( double *f , const void *data , const double *fparams )
{
  double *y = ((struct data *)data)->y;
  const double A = fparams[0] , lambda = fparams[1] , b = fparams[2] ;
  size_t i ,n = ((struct data *)data)->n ;
  for (i = 0; i < n; i++) {
    f[i] = A * exp (-lambda * i) + b - y[i] ;
  }
  return ;
}

// derivatives
void
exp_df( double **df , const void *data , const double *fparams )
{
  const double A = fparams[0] , lambda = fparams[1] ;
  size_t i , n = ((struct data *)data)->n ;
  for( i = 0 ; i < n ; i++ ) {
    const double t = i;
    const double e = exp(-lambda * t);
    df[0][i] = e ;
    df[1][i] = -t * A * e ;
    df[2][i] = 1.0 ;
  }
  return ;
}

// second derivatives
void
exp_d2f( double **d2f , const void *data , const double *fparams )
{
  const double A = fparams[0] , lambda = fparams[1] ;
  size_t i , n = ((struct data *)data)->n ;
  for( i = 0 ; i < n ; i++ ) {
    const double t = i;
    const double e = exp(-lambda * t);
    d2f[0][i] = 0.0    ; d2f[1][i] = -t * e  ; d2f[2][i] = 0.0 ;
    d2f[3][i] = -t * e ; d2f[4][i] = t*t*A*e ; d2f[5][i] = 0.0 ;
    d2f[6][i] = 0.0    ; d2f[7][i] = 0.0     ; d2f[8][i] = 0.0 ;
  }
  return ;
}

void
exp_guesses( double *fparams )
{
  if( fparams[0] == UNINIT_FLAG && fparams[1] == UNINIT_FLAG &&
      fparams[2] == UNINIT_FLAG ) {
    fparams[0] = 7.0 ; fparams[1] = 0.2 ; fparams[2] = 1.1 ;
  }
  return ;
}

void
exp_priors( double *priors , double *err_priors )
{
  /*
  priors[0] = 5.0 ; err_priors[0] = 0.1 ;
  priors[1] = 0.1 ; err_priors[1] = 0.005 ;
  priors[2] = 1.0 ; err_priors[2] = 0.05 ;
  */
}
