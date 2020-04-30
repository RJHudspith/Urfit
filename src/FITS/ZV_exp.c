/**
   exponential Model f = ( ZV * ( 1 - exp(-lambda * i) )
 */
#include "gens.h"

double
fZV_exp( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[0] * ( 1 - exp( -fparams[1] * X.X ) ) ;
}

void
ZV_exp_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i ;
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] , DATA->N, DATA->M } ;
    f[i] = fparams[ DATA -> map[i].p[0] ] * ( 1 - exp( -fparams[ DATA -> map[i].p[1] ] * X.X ) ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
ZV_exp_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = ( const struct data* )data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double e = exp(-fparams[DATA -> map[i].p[1]] * t) ;
    df[DATA -> map[i].p[0]][i] = 1 - e ;
    df[DATA -> map[i].p[1]][i] = +fparams[DATA -> map[i].p[0]]*t*e ;
  }
  return ;
}

// second derivatives
void
ZV_exp_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
ZV_exp_guesses( double *fparams ,
		const struct data_info Data ,
		const struct fit_info Fit )
{
  // just set them to one
  size_t i ;
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    fparams[0] = 1.0 ;
  }
  
  return ;
}
