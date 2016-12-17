/**
   PADE Model f = a + ( b*x + ... + b_n*x^n) / ( 1 + c*x + ... + c_mx^m )
 */
#include "gens.h"

// set the n and m
static size_t n = 1 , m = 1 ;

void 
pade_set_nm( const size_t new_n ,
	     const size_t new_m )
{
  if( new_n > 0 && new_m > 0 ) {
    n = new_n ; m = new_m ;
  }
}

void
pade_get_nm( size_t *new_n ,
	     size_t *new_m )
{
  *new_n = n ; *new_m = m ;
}


double
fpade( const struct x_desc X , const double *fparams , const size_t Npars )
{
  register double numerator = X.X * fparams[ n ] ;
  size_t i ;
  for( i = n-1 ; i > 0 ; i-- ) {
    numerator = X.X * ( fparams[i] + numerator ) ;
  }
  register double denominator = X.X * fparams[ n + m ] ;
  for( i = n + m - 1 ; i > n ; i-- ) {
    denominator = X.X * ( fparams[i] + denominator ) ;
  }
  return fparams[0] + ( numerator ) / ( 1 + denominator ) ;
}

void
pade_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT } ;
    f[i] = fpade( X , fparams , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
pade_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double X = DATA -> x[i] ;
    // precompute numerator
    double numerator = X * fparams[ n ] ;
    for( j = n - 1 ; j > 0 ; j-- ) {
      numerator = X * ( fparams[j] + numerator ) ;
    }
    // precompute denominator
    register double denominator = X * fparams[ n + m ] ;
    for( j = n + m - 1 ; j > n ; j-- ) {
      denominator = X * ( fparams[j] + denominator ) ;
    }
    denominator += 1.0 ;

    // this factor is pretty common
    const double factor = -( numerator ) / ( denominator * denominator ) ;

    df[0][i] = 1.0 ;
    // numerator derivatives
    for( j = 1 ; j <= n ; j++ ) {
      df[j][i] = pow( X , j ) / denominator ;
    }
    // denominator  derivatives
    for( j = 0 ; j < m ; j++ ) {
      df[j+n+1][i] = pow( X , j+1 ) * factor ;
    }
  }
  return ;
}

// second derivatives
void
pade_d2f( double **d2f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    for( j = 0 ; j < (n+m+1)*(n+m+1) ; j++ ) {
      // these derivatives are hard
      d2f[j][i] = 0.0 ;
    }
  }
  return ;
}

void
pade_guesses( double *fparams )
{
  return ;
}

void
pade_priors( double *priors , double *err_priors )
{
  return ;
}
