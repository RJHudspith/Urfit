/**
   multi-exponential Model plus constant 
 */
#include "gens.h"

double
fexp_plusc( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  register double sum = fparams[0] ;
  for( i = 0 ; i < 2*X.N ; i+=2 ) {
    sum += fparams[i+1] * exp( -fparams[i+2] * X.X ) ;
  }
  return sum ;
}

void
exp_plusc_f( double *f , const void *data , const double *fparams )
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
    f[i] = fexp_plusc( X , p , DATA -> N*2+1 ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
exp_plusc_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    // constant
    df[ DATA-> map[i].p[0] ][i] = 1 ;
    // exponentials
    for( j = 0 ; j < DATA -> N*2 ; j+=2 ) {
      const double t = DATA -> x[i] ;
      const double e = exp( -fparams[ DATA-> map[i].p[j+2] ] * t ) ;
      df[ DATA-> map[i].p[j+1] ][i] = e ;
      df[ DATA-> map[i].p[j+2] ][i] = -t * fparams[ DATA-> map[i].p[j+1] ] * e ;
    }
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
exp_plusc_guesses( double *fparams , const size_t Nlogic )
{
  size_t i , all_flagged = 0 ;
  for( i = 0 ; i < Nlogic ; i++ ) {
    if( fparams[i] == UNINIT_FLAG ) {
      all_flagged++ ;
    }
  }

  // perform a guess, otherwise assume someone has set them
  if( all_flagged == Nlogic ) {
    fparams[0] = 1.0 ; fparams[1] = 200.0 ; 
    fparams[2] = 1.0 ; fparams[3] = 1.0 ;
  }
  
  return ;
}
