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
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT } ;
    //f[i] = fexp( X , fparams , DATA -> Npars ) - DATA -> y[i] ;
    f[i] = fparams[ DATA-> map[i].p[1] ] *
      exp( -fparams[ DATA-> map[i].p[0] ] * X.X ) - DATA -> y[i] ;
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
    const double e = exp(-fparams[ DATA-> map[i].p[0] ] * t ) ;
    df[ DATA-> map[i].p[0] ][i] = -t * fparams[ DATA-> map[i].p[1] ] * e ;
    df[ DATA-> map[i].p[1] ][i] = e ;
  }
  return ;
}

// second derivatives
void
exp_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
exp_guesses( double *fparams )
{
  if( fparams[0] == UNINIT_FLAG && fparams[1] == UNINIT_FLAG ) {
    fparams[0] = 1.0 ; fparams[1] = 200.0 ; 
    fparams[2] = 1.0 ; fparams[3] = 1.0 ; 
  }
  return ;
}

void
exp_priors( double *priors , double *err_priors )
{
  return ;
}
