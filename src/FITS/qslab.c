/**
   @file qslab.c
   @brief fit the slab formula Y = \chi X + C * ( 1 - e^{m t_ref}) * ( 1 - e^{m( L_T - t_ref ) ) 
 */
#include "gens.h"
#include "Nder.h"

double
fqslab( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[0] * X.X +
    fparams[1] * ( 1.0 - exp( -fparams[2] * X.X ) - exp( -fparams[2] * ( X.LT - X.X ) ) + exp( -fparams[2] * X.LT ) ) ;
}

void
qslab_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0; i < DATA -> n ; i++ ) {
    double p[ DATA -> Npars ] ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;
    f[i] = fqslab( X , p , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
qslab_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double x = DATA -> x[i] ;
    const double m = fparams[ DATA -> map[i].p[2] ] ;
    const double C = fparams[ DATA -> map[i].p[1] ] ;
    const double LT = DATA -> LT[i] ;
    
    df[ DATA -> map[i].p[0] ][ i ] = x ;
    df[ DATA -> map[i].p[1] ][ i ] = \
      ( 1.0 - exp( -m * x ) - exp( -m * ( LT - x ) ) + exp( -m * LT ) ) ;
    df[ DATA -> map[i].p[2] ][ i ] = \
      C * ( x * exp( -m * x ) + ( LT - x ) * exp( -m * ( LT - x ) ) - LT * exp( -m * LT ) ) ;
  }
  return ;
}

// second derivatives
void
qslab_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// guesses using the data, we don't need to be too accurate
void
qslab_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{ 
  fparams[0] = Data.y[ Data.Ndata[0] - 1 ].avg / Data.LT[ Data.Ndata[0] - 1 ] ;
  fparams[1] = -1.385E-04 ;
  fparams[2] = -0.167 ;
}
