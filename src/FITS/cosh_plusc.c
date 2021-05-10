/**
   multi-exponential Model plus constant 
 */
#include "gens.h"

double
fcosh_plusc( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  register double sum = fparams[0] ;
  for( i = 0 ; i < 2*X.N ; i+=2 ) {
    sum +=
      fparams[i+1] * ( exp( -fparams[i+2] * X.X ) + exp( -fparams[i+2] * ( X.LT/2. - X.X ) ) ) ;
  }
  return sum ;
}

void
cosh_plusc_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i , j ;
  for (i = 0; i < DATA -> n ; i++) {
    double p[ DATA -> N*2+1 ] ;
    for( j = 0 ; j < DATA -> N*2+1 ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fcosh_plusc( X , p , DATA -> N*2+1 ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
cosh_plusc_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    // constant
    df[ DATA-> map[i].p[0] ][i] = 1 ;
    // exponentials
    for( j = 0 ; j < DATA -> N*2 ; j+=2 ) {
      const double t  = DATA -> x[i] ;
      const double e1 = exp( -fparams[ DATA-> map[i].p[j+2] ] * t ) ;
      const double e2 = exp( -fparams[ DATA-> map[i].p[j+2] ] * ( DATA -> LT[i]/2. - t ) ) ;
      df[ DATA-> map[i].p[j+1] ][i] = e1+e2 ;
      df[ DATA-> map[i].p[j+2] ][i] =
	fparams[ DATA-> map[i].p[j+1] ] * (-t * e1 -( DATA -> LT[i]/2. - t )*e2 ) ;
    }
  }
  return ;
}

// second derivatives
void
cosh_plusc_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
cosh_plusc_guesses( double *fparams , const size_t Nlogic )
{
  // do nothing
  return ;
}
