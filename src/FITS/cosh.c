/**
   multiple cosh model

   f = sum_n A < exp( -lambda * x[i] ) + exp( -lambda * ( Lt - x[i] ) ) >
 */
#include "gens.h"

double
fcosh( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  register double sum = 0.0 ;
  for( i = 0 ; i < 2 * X.N ; i+=2 ) {
    sum += fparams[i] * ( exp( -fparams[i+1] * X.X ) +
			  exp( -fparams[i+1] * ( X.LT - X.X ) ) ) ;
  }
  return sum ;
}

void
cosh_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0; i < DATA -> n ; i++ ) {
    double p[ DATA -> N * 2 ] ;
    for( j = 0 ; j < DATA -> N * 2 ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;
    f[i] = fcosh( X , p , DATA -> N * 2 ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
cosh_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    for( j = 0 ; j < 2*DATA -> N ; j+= 2 ) {
      const double t = DATA -> x[i] ;
      const double fwd = exp( -fparams[ DATA -> map[i].p[j+1] ] * t ) ;
      const double bwd = exp( -fparams[ DATA -> map[i].p[j+1] ] * ( DATA -> LT[i] - t ) ) ;
      df[ DATA -> map[i].p[j+0] ][i] = fwd + bwd ;
      df[ DATA -> map[i].p[j+1] ][i] = -fparams[ DATA -> map[i].p[j] ] *
	( t * fwd + ( DATA -> LT[i] - t ) * bwd ) ;
    }
  }
  return ;
}

// second derivatives
void
cosh_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}
