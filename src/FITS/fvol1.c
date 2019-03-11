/**
   @file fvol1.c
   @brief finite volume fit m_pi + A * e^{-m_pi X}
 */
#include "gens.h"

#include "Nder.h"

size_t L[3] = { 32 , 48 , 64 } ;

double
ffvol1( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[0] + fparams[1] * exp( -fparams[0] * X.X ) ;
}

void
fvol1_f( double *f , const void *data , const double *fparams )
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
    f[i] = ffvol1( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvol1_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double X = DATA -> x[i] ;
    const double expfac = exp( -fparams[ DATA -> map[i].p[0] ] * X ) ;
    df[0][i] = 1.0 - X * fparams[ DATA -> map[i].p[1] ] * expfac ;
    df[1][i] = expfac ;
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvol1_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvol1_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 0.05 ;
  fparams[1] = 1 ;
  //fparams[2] = -400 ;
  return ;
}
