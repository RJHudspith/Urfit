/**
   exponential Model f = ( A * exp(-lambda * i) - y_i )
 */
#include "gens.h"

double
fcosh( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[1] * ( exp( -fparams[0] * X.X ) + 
			exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
}

void
cosh_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT } ;
    f[i] = fcosh( X , fparams , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
cosh_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double fwd = exp(-fparams[0] * t ) ;
    const double bwd = exp(-fparams[0] * ( DATA -> LT - t ) ) ;
    df[0][i] = -fparams[1] * ( t * fwd + ( DATA -> LT - t ) * bwd ) ;
    df[1][i] = fwd + bwd ;
  }
  return ;
}

// second derivatives
void
cosh_d2f( double **d2f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double fwd = exp(-fparams[0] * t ) ;
    const double bwd = exp(-fparams[0] * ( DATA -> LT - t ) ) ;
    d2f[0][i] = -fparams[1] * ( t * t * fwd + 
				( DATA -> LT - t ) * ( DATA -> LT - t ) * bwd ) ;
    d2f[1][i] = -( t * fwd + ( DATA -> LT - t ) * bwd ) ;
    d2f[2][i] = d2f[1][i] ;
    d2f[3][i] = 0.0 ;
  }
  return ;
}

void
cosh_guesses( double *fparams )
{
  if( fparams[0] == UNINIT_FLAG && fparams[1] == UNINIT_FLAG ) {
    fparams[0] = 0.2 ; fparams[1] = 7.0 ;
  }
  return ;
}

void
cosh_priors( double *priors , double *err_priors )
{
  return ;
}
