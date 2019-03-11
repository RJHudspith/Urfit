/**
   single tanh model f = p[0]*(exp[-p[1]*t] - exp(-p[1]*(T-t)))/(exp[-p[1]*t] + exp(-p[1]*(T-t)))
 */
#include "gens.h"

#include "Nder.h"

double
ftanh( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double fwd = exp( -fparams[1] * X.X ) ;
  const double bwd = exp( -fparams[1] * ( X.LT - X.X ) ) ;
  return fparams[0] * ( fwd - bwd ) / ( fwd + bwd ) ;
}

void
tanh_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for (i = 0; i < DATA -> n ; i++) {
    double p[ 2 ] ;
    for( j = 0 ; j < 2 ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;
    f[i] = ftanh( X , p , DATA -> N * 2 ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
tanh_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( ftanh , X ,
					   DATA -> map[i].bnd ,
					   fparams , j , DATA -> Npars ) ;
    }
  }
  return ;
}

// second derivatives
void
tanh_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

