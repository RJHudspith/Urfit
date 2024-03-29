/**
   @file poles.c
   @brief fit some poles
 */
#include "gens.h"

#include "Nder.h"

double
fpoles( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  register double sum = 0.0 ;
  for( i = 0 ; i < X.N ; i++ ) {
    sum += fparams[ i ] * pow( X.X , -(int)X.N + (int)i ) ; 
  }
  for( i = 0 ; i < X.M+1 ; i++ ) {
    sum += fparams[ X.N + i ] * pow( X.X , i ) ;
  }
  return sum ;
  //return fparams[0]/(X.X) + fparams[1] + fparams[2]*X.X ;
}

void
poles_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> Npars ] ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fpoles( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
poles_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    for( j = 0 ; j < X.N ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = fparams[ j ] * pow( X.X , -(int)X.N + (int)j ) ; 
    }
    for( j = 0 ; j < X.M+1 ; j++ ) {
      df[ DATA -> map[i].p[j+X.N] ][i] = fparams[ j+X.N ] * pow( X.X , (int)X.N + (int)j ) ; 
    }
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
poles_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// compute the linearised matrix for the GLS
void
poles_linmat( double **U ,
	      const void *data ,
	      const size_t N ,
	      const size_t M ,
	      const size_t Nlogic )
{
  struct data *Data = (struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < Data -> n ; i++ ) {
    // set matrix to zero
    for( j = 0 ; j < Nlogic ; j++ ) {
      U[i][j] = 0.0 ;
    }
    const double x0 = Data -> x[i] ;
    for( j = 0 ; j < N ; j++ ) {
      U[i][ Data -> map[i].p[j] ] = pow( x0 , -(int)N+(int)j ) ;
    }
    for( j = 0 ; j < M+1 ; j++ ) {
      U[i][ Data -> map[i].p[j + N] ] = pow( x0 , (int)j ) ;
    }
  }
  return ;
}

void
poles_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  // do something here
  size_t i ;
  for( i = 0 ; i < Fit.N + Fit.M + 1 ; i++ ) {
    fparams[i] = 1.0 ;
  }
}
