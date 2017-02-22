/**
   exponential Model f = ( A * exp(-lambda * i) - y_i )
 */
#include "gens.h"

double
fsinh( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[0] * ( exp( -fparams[1] * X.X ) - 
			exp( -fparams[1] * ( X.LT - X.X ) ) ) ;
}

void
sinh_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] } ;
    f[i] = fparams[ DATA -> map[i].p[0] ] *
      ( exp( -fparams[ DATA -> map[i].p[1] ] * X.X ) - 
	exp( -fparams[ DATA -> map[i].p[1] ] * ( X.LT - X.X ) ) )
      - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
sinh_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double fwd = exp(-fparams[ DATA -> map[i].p[1] ] * t ) ;
    const double bwd = exp(-fparams[ DATA -> map[i].p[1] ] * ( DATA -> LT[i] - t ) ) ;
    df[  DATA -> map[i].p[1] ][i] = -fparams[ DATA -> map[i].p[0] ] *
      ( t * fwd - ( DATA -> LT[i] - t ) * bwd ) ;
    df[  DATA -> map[i].p[0] ][i] = fwd - bwd ;
  }
  return ;
}

// second derivatives
void
sinh_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

