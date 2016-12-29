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
    f[i] = //fexp_plusc( X , fparams , DATA -> Npars )
      fparams[ DATA -> map[i].p[1] ] *
      exp( -fparams[ DATA -> map[i].p[0] ] * X.X ) +
      fparams[ DATA -> map[i].p[2] ]
      - DATA -> y[i] ;
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
    const double e = exp(-fparams[DATA -> map[i].p[0]] * t) ;
    df[DATA -> map[i].p[0]][i] = -t * fparams[DATA -> map[i].p[1]] * e ;
    df[DATA -> map[i].p[1]][i] = e ;
    df[DATA -> map[i].p[2]][i] = 1.0 ;
  }
  return ;
}

// second derivatives
void
exp_plusc_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
exp_plusc_guesses( double *fparams )
{
  if( fparams[0] == UNINIT_FLAG || fparams[1] == UNINIT_FLAG ||
      fparams[2] == UNINIT_FLAG ) {
    fparams[0] = 0.12 ; fparams[1] = 5.5 ; fparams[2] = 10 ;
  }
  return ;
}

void
exp_plusc_priors( double *priors , double *err_priors )
{
}
