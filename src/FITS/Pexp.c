/**
   HadSpec principle correlator fit
 */
#include "gens.h"

double
fPexp( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double e = exp( -fparams[0] * X.X ) ;
  double sum = 0.0 ;
  size_t i ;
  for( i = 0 ; i < Npars ; i++ ) {
    sum += fparams[1+2*i] * ( exp( -fparams[2+2*i] * X.X ) - e ) ;
  }
  return e + sum ;
}

void
Pexp_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i , j ;
  for (i = 0; i < DATA -> n ; i++) {
    double p[ 1+2*DATA->N ] ;
    for( j = 0 ; j < 1+2*DATA->N ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = ( fPexp( X , p , DATA->N ) - DATA -> y[i] ) ;
  }
  return ;
}

// derivatives
void
Pexp_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e1 = exp(-fparams[DATA -> map[i].p[0]] * t) ;
    size_t j ;
    for( j = 0 ; j < DATA->N ; j++ ) {
      const double A = fparams[DATA -> map[i].p[1+2*j]] ;
      const double e2 = exp(-fparams[DATA -> map[i].p[2+2*j]] * t) ;
      df[DATA -> map[i].p[0]][i] = -t * ( 1 + A ) * e1 ;
      df[DATA -> map[i].p[1]][i] = -e1 + e2 ;
      df[DATA -> map[i].p[2]][i] = t * A * fparams[DATA -> map[i].p[2+2*j]] * e2 ;
    }
  }
  return ;
}

// second derivatives
void
Pexp_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
Pexp_guesses( double *fparams , const size_t Nlogic )
{
  size_t i , all_flagged = 0 ;
  for( i = 0 ; i < Nlogic ; i++ ) {
    if( fparams[i] == UNINIT_FLAG ) {
      all_flagged++ ;
    }
  }

  // perform a guess, otherwise assume someone has set them
  if( all_flagged == Nlogic ) {
    fparams[0] = 1.0 ; fparams[1] = 1.0 ; fparams[2] = 1.0 ; 
  }
  
  return ;
}
