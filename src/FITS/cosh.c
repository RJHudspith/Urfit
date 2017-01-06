/**
   exponential Model f = ( A * exp(-lambda * i) - y_i )
 */
#include "gens.h"

double
fcosh( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[1] * ( exp( -fparams[0] * X.X ) + 
			exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
}

void
cosh_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT } ;

    f[i] = fparams[ DATA -> map[i].p[1] ] *
      ( exp( -fparams[ DATA -> map[i].p[0] ] * X.X ) + 
	exp( -fparams[ DATA -> map[i].p[0] ] * ( X.LT - X.X ) ) )
      - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
cosh_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double fwd = exp( -fparams[ DATA -> map[i].p[0] ] * t ) ;
    const double bwd = exp( -fparams[ DATA -> map[i].p[0] ] * ( DATA -> LT - t ) ) ;
    df[ DATA -> map[i].p[0] ][i] = -fparams[ DATA -> map[i].p[1] ] *
      ( t * fwd + ( DATA -> LT - t ) * bwd ) ;
    df[ DATA -> map[i].p[1] ][i] = fwd + bwd ;
  }
  return ;
}

// second derivatives
void
cosh_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
cosh_guesses( double *fparams , const size_t Nlogic )
{
  size_t i , all_flagged = 0 ;
  for( i = 0 ; i < Nlogic ; i++ ) {
    if( fparams[i] == UNINIT_FLAG ) {
      all_flagged++ ;
    }
  }

  // perform a guess, otherwise assume someone has set them
  if( all_flagged == Nlogic ) {
    for( i = 0 ; i < Nlogic ; i++ ) {
      fparams[i] = 1 + i ;
    }
  }
  
  return ;
}

void
cosh_priors( double *priors , double *err_priors )
{
  return ;
}
