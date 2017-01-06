/**
   exponential Model f = ( A * exp(-lambda * i) + b - y_i )
 */
#include "gens.h"

double
fexp_plusc( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[1] * exp( -fparams[0] * X.X ) + fparams[2] ;
}

void
exp_plusc_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i ;
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA-> LT } ;
    f[i] = //fexp_plusc( X , fparams , DATA -> Npars )
      fparams[ DATA -> map[i].p[1] ] *
      exp( -fparams[ DATA -> map[i].p[0] ] * X.X ) +
      fparams[ DATA -> map[i].p[2] ]
      - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
exp_plusc_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e = exp(-fparams[DATA -> map[i].p[0]] * t) ;
    df[DATA -> map[i].p[0]][i] = -t * fparams[DATA -> map[i].p[1]] * e ;
    df[DATA -> map[i].p[1]][i] = e ;
    df[DATA -> map[i].p[2]][i] = 1.0 ;
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

void
exp_plusc_priors( double *priors , double *err_priors )
{
}
