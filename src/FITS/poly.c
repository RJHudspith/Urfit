/**
   exponential Model f = ( A * exp(-lambda * i) - y_i )
 */
#include "gens.h"

static size_t n = 0 ;

void
poly_set_n( const size_t new_n )
{
  n = new_n ;
  return ;
}

double
fpoly( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  if( Npars < 2 ) return fparams[0] ;
  register double poly = X.X * fparams[ Npars-1 ] ;
  for( i = Npars-2 ; i > 0 ; i-- ) {
    poly = X.X * ( fparams[i] + poly ) ;
  }
  return poly + fparams[0] ;
}

void
poly_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , p ; 
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT } ;
    double par[ DATA -> Npars ] ;
    for( p = 0 ; p < DATA -> Npars ; p++ ) {
      par[ p ] = fparams[ DATA -> map[i].p[p] ] ;
    }
    f[i] = fpoly( X , par , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
poly_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double xloc = 1.0 ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = xloc ;
      xloc *= DATA->x[i] ;
    }
  }
  return ;
}

// second derivatives
void
poly_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
poly_guesses( double *fparams , const size_t Nlogic )
{
  // perform a good guess
  size_t i , all_flagged = 0 ;
  for( i = 0 ; i < Nlogic ; i++ ) {
    if( fparams[i] == UNINIT_FLAG ) {
      all_flagged++ ;
    }
  }

  // perform a guess, otherwise assume someone has set them
  if( all_flagged == Nlogic ) {
    for( i = 0 ; i < Nlogic ; i++ ) {
      fparams[i] = i+1 ;
    }
  }
  return ;
}

void
poly_priors( double *priors , double *err_priors )
{
  return ;
}