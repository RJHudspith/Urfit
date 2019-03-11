/**
   Fits f[0]*f[1]*exp^{-f[2]t}
 */
#include "gens.h"

double
fHALexp( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  register double sum = 0.0 ;
  for( i = 0 ; i < 3 * X.N ; i+=3 ) {
    sum += fparams[i] * fparams[i+1] * exp( -fparams[i+2] * X.X ) ;
  }
  return sum ;
}

void
HALexp_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> N * 3 ] ;
    for( j = 0 ; j < DATA -> N * 3 ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fHALexp( X , p , DATA -> N * 3 ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
HALexp_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {    
    for( j = 0 ; j < 3*DATA -> N ; j+=3 ) {
      const double t = DATA -> x[i] ;
      const double e = exp( -fparams[ DATA-> map[i].p[j+2] ] * t ) ;
      df[ DATA-> map[i].p[j+0] ][i] = fparams[ DATA-> map[i].p[j+1] ] * e ;
      df[ DATA-> map[i].p[j+1] ][i] = fparams[ DATA-> map[i].p[j+0] ] * e ;
      df[ DATA-> map[i].p[j+2] ][i] = -t * fparams[ DATA-> map[i].p[j+0] ] * fparams[ DATA-> map[i].p[j+1] ] * e ;
    }
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
HALexp_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// guesses using the data, we don't need to be too accurate
void
HALexp_guesses( double *fparams ,
		const struct data_info Data ,
		const struct fit_info Fit )
{
  return ;
}
