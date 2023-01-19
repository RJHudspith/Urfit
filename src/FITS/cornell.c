/**
   @file cornell.c
   @brief cornell fit for the SU(2) project
 */
#include "gens.h"

#ifdef MDEP
// this is the bare m-m^\chiral
static const double mbar[ 8 ] = { 0.056 , // 0.845
				  0.046 , // 0.855
				  0.036 , // 0.865
				  0.026 , // 0.875
				  0.021 , // 0.880
				  0.016 , // 0.885
				  0.011 , // 0.890
				  0.001 , // 0.900
} ;
#endif

double
fcornell( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef MDEP
  return
      fparams[0] * ( 1 + mbar[ Npars ] * fparams[1] ) * X.X 
    + fparams[2] * ( 1 + mbar[ Npars ] * fparams[3] )
    + fparams[4] * ( 1 + mbar[ Npars ] * fparams[5] ) / X.X ;
#else
  return fparams[0] * X.X  + fparams[1] + fparams[2] / X.X ;
#endif
}

void
cornell_f( double *f , const void *data , const double *fparams )
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
    f[i] = fcornell( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

void
cornell_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double r = DATA -> x[i] ;
#ifdef MDEP
    const double m = mbar[ DATA -> map[i].bnd ];
    // derivatives wrt the various fit params
    df[0][i] = ( 1 + m * fparams[1] ) * r ;
    df[1][i] = ( fparams[0] * m ) * r ;
    df[2][i] = ( 1 + m * fparams[3] ) ;
    df[3][i] = ( fparams[2] * m ) ;
    df[4][i] = ( 1 + m * fparams[5] ) / r ;
    df[5][i] = ( fparams[4] * m ) / r ;
#else
    df[0][i] = r ;
    df[1][i] = 1 ;
    df[2][i] = 1./r ;
#endif
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
cornell_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
cornell_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
#ifdef MDEP
  fparams[0] = 0.12 ;
  fparams[1] = 0.1 ;
  fparams[2] = 1.0 ;
  fparams[3] = 1.0 ;
  fparams[4] = 0.1 ;
  fparams[5] = 0.01 ;
#else
  fparams[0] = -1 ;
  fparams[1] = 0.5 ;
  fparams[2] = 0.1 ;
#endif
}
