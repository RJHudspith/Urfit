/**
   Asymmetric multi cosh model for instance a baryon looks like this

   f = \sum_n A_n exp( -lambda_n * x[i] ) + \sum_m A_m exp( -lambda_m * ( Lt - x[i] ) )
 */
#include "gens.h"

double
fcosh_asymm( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  register double sum = 0.0 ;
  for( i = 0 ; i < 2*X.N ; i+=2 ) {
    sum += fparams[i] * exp( -fparams[i+1] * X.X ) ;
  }
  for( i = 2*X.N ; i < 2*(X.N+X.M) ; i+=2 ) {
    sum += fparams[i]* exp( -fparams[i+1] * ( X.LT - X.X ) ) ;
  }
  return sum ;
}

void
cosh_asymm_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0; i < DATA -> n ; i++ ) {
    double p[ 2*(DATA->N+DATA->M) ] ;
    for( j = 0 ; j < 2*(DATA -> N + DATA->M) ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] , DATA -> N , DATA -> M } ;
    f[i] = fcosh_asymm( X , p , 2*(DATA -> N + DATA->M) ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
cosh_asymm_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    for( j = 0 ; j < 2*DATA -> N ; j+= 2 ) {
      const double fwd = exp( -fparams[ DATA -> map[i].p[j+1] ] * t ) ;
      df[ DATA -> map[i].p[j+0] ][i] = fwd ;
      df[ DATA -> map[i].p[j+1] ][i] = -fparams[ DATA -> map[i].p[j+0] ] *( t * fwd ) ;
    }
    for( j = 2*DATA->N ; j < 2*(DATA -> N+DATA->M) ; j+= 2 ) {
      const double bwd = exp( -fparams[ DATA -> map[i].p[j+1] ] * ( DATA -> LT[i] - t ) ) ;
      df[ DATA -> map[i].p[j+0] ][i] = bwd ;
      df[ DATA -> map[i].p[j+1] ][i] = -fparams[ DATA -> map[i].p[j+0] ] * ( DATA -> LT[i] - t ) * bwd ;
    }
  }
  return ;
}

// second derivatives
void
cosh_asymm_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}
