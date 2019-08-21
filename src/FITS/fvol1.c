/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

size_t L[3] = { 48 , 48 , 48 } ;

double
ffvol1( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double MPIL = sqrt(X.X)*48 ;
  return fparams[0] + fparams[1]*(X.X) + fparams[2]*exp( -MPIL ) ; ///sqrt(MPIL) ;
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
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( ffvol1 , X ,
					   DATA -> map[i].bnd ,
					   fparams , j , DATA -> Npars ) ;
    }
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
  fparams[2] = -400 ;
  return ;
}
