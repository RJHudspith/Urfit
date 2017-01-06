/**
   exponential Model f = ( A * exp(-lambda * i) - y_i )
 */
#include "gens.h"

double
fexp( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[1] * exp( -fparams[0] * X.X ) ;
}

void
exp_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT } ;
    f[i] = fparams[ DATA-> map[i].p[1] ] *
      exp( -fparams[ DATA-> map[i].p[0] ] * X.X ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
exp_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e = exp(-fparams[ DATA-> map[i].p[0] ] * t ) ;
    df[ DATA-> map[i].p[0] ][i] = -t * fparams[ DATA-> map[i].p[1] ] * e ;
    df[ DATA-> map[i].p[1] ][i] = e ;
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
exp_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
exp_guesses( double *fparams , const size_t Nlogic )
{
  bool flag = false ;
  size_t i ;
  for( i = 0 ; i < Nlogic ; i++ ) {
    if( fparams[i] != UNINIT_FLAG ) {
      flag = true ;
      continue ;
    }
  }

  // perform a guess, otherwise assume someone has set them
  if( flag == false ) {
    for( i = 0 ; i < Nlogic ; i++ ) {
      fparams[i] = i + 1 ;
    }
  }

  return ;
}

void
exp_priors( double *priors , double *err_priors )
{
  return ;
}
