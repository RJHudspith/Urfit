/**
   @file Qcorr_Bessel.c
   @brief computes the topological correlator mass using

   y = A/r^1.5*(1+3/(8*m*r))exp(-mr)
 */
#include "gens.h"

double
fQcorr_bessel( const struct x_desc X , const double *fparams , const size_t Npars )
{
  double sum = 0 ;
  for( int n = 0 ; n < 2*X.N ; n+=2 ) { 
    sum += fparams[n]*(1+3/(8*fparams[n+1]*X.X))*exp(-fparams[n+1]*X.X);
  }
  return sum/(X.X*sqrt(X.X)) ;
}

void
Qcorr_bessel_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> N * 2 ] ;
    for( j = 0 ; j < DATA -> N * 2 ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fQcorr_bessel( X , p , DATA -> N * 2 ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
Qcorr_bessel_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    const double rfac = (1./(X.X*sqrt(X.X))) ;
    for( j = 0 ; j < 2*DATA->N ; j+=2 ) {

      const double A = fparams[DATA -> map[ i ].p[j+0]] ;
      const double m = fparams[DATA -> map[ i ].p[j+1]] ;
      const double E = exp( -m*X.X ) ;
      const double P = (1+3/(8*m*X.X)) ;
      
      df[ DATA -> map[ i ].p[j+0] ][i] = P*E*rfac ;
      df[ DATA -> map[ i ].p[j+1] ][i] = -A*E*(X.X+3/(8*m)+3/(8*m*m*X.X))*rfac ;
      
    }
  }
  return ;
}

void
Qcorr_bessel_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
Qcorr_bessel_guesses( double *fparams ,
	     const struct data_info Data ,
	     const struct fit_info Fit )
{
  return ;
}
