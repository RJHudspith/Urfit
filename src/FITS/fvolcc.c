/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

//#define V3

// TODO
/*
(A+B a2 MetaC^2)/(
1/4 (-C + MetaC)^2 + (D+
   E a2 MetaC^2)^2)
 */

double
ffvolcc( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double asq = fparams[ 4 + Npars ] ;

#ifdef V2
  return
    fparams[0] * X.X * ( 1. + fparams[1] * asq * X.X + fparams[2] * asq * asq * X.X * X.X )
    + fparams[3]*asq ;
#elif (defined V3)
  return fparams[0]*asq
    +fparams[1]*X.X
    +fparams[2]*asq*asq
    +fparams[3]*asq*asq*X.X
    ;
#elif (defined V4)
  return fparams[0]*asq
    +fparams[1]*X.X
    +fparams[2]*asq*asq
    +fparams[3]*asq*X.X
    ;
#else
  return fparams[0]
    +fparams[1]*X.X
    +fparams[2]*asq
    +fparams[3]*asq*X.X
    ;
#endif
}

void
fvolcc_f( double *f , const void *data , const double *fparams )
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
    f[i] = ffvolcc( X , p , DATA->map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvolcc_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    const double asq = fparams[ 4 + DATA->map[i].bnd ] ;
    
#ifdef V2
    df[ DATA -> map[i].p[0] ][i] = X.X * ( 1 + fparams[1] * asq * X.X + fparams[2] * asq * asq * X.X * X.X ) ;
    df[ DATA -> map[i].p[1] ][i] = fparams[0]*X.X * asq * X.X ;
    df[ DATA -> map[i].p[2] ][i] = fparams[0]*X.X * asq * asq * X.X * X.X ;
    df[ DATA -> map[i].p[3] ][i] = asq ;
    df[ 4 + DATA->map[i].bnd ][ i ] =
      fparams[0] * X.X*X.X*( fparams[1] + fparams[2]*X.X*2*asq ) + fparams[3] ;
#elif (defined V3)
    df[ DATA -> map[i].p[0] ][i] = asq ;
    df[ DATA -> map[i].p[1] ][i] = X.X ;
    df[ DATA -> map[i].p[2] ][i] = asq*asq ;
    df[ DATA -> map[i].p[3] ][i] = asq*asq*X.X ;
    df[ 4 + DATA->map[i].bnd ][ i ] =
      fparams[0] + 2*asq*fparams[2] + 2*asq*fparams[3]*X.X ;
#elif (defined V4)
    df[ DATA -> map[i].p[0] ][i] = asq ;
    df[ DATA -> map[i].p[1] ][i] = X.X ;
    df[ DATA -> map[i].p[2] ][i] = asq*asq ;
    df[ DATA -> map[i].p[3] ][i] = asq*X.X ;
    df[ 4 + DATA->map[i].bnd ][ i ] =
      fparams[0] + 2*asq*fparams[2] + fparams[3]*X.X ; 
#else
    df[ DATA -> map[i].p[0] ][i] = 1.0 ;
    df[ DATA -> map[i].p[1] ][i] = X.X ;
    df[ DATA -> map[i].p[2] ][i] = asq ;
    df[ DATA -> map[i].p[3] ][i] = asq*X.X ;
    df[ 4 + DATA->map[i].bnd ][ i ] =
      fparams[2] + fparams[3]*X.X ;
#endif
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvolcc_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvolcc_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 0.0 ;
  fparams[1] = 1 ;
  return ;
}
