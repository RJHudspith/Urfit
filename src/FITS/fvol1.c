/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

//#define L3_EXTRAP
//#define MPIL

#define V1

#ifdef MKL
// U103, H101, U102 , H102, H105, N101, C101
static const double MPIL[13] = {
  4.3269264 , 5.8321216 ,
  4.6145688 , 6.100112 ,
  6.436992 , 9.6920208 ,
  9.8667504 ,
  4.5 , 5.5 , 6.5 , 1E8 // idea is to have the infinite volume guy here
} ;

#else
// U103, H101, U102 , H102, U101, H105, N101, C101
static const double MPIL[13] = {
  4.3269264 , 5.8321216 ,
  3.7376688 , 4.9293248,
  2.9499072 , 3.91064 , 5.887488 ,
  4.575890879999999 ,
  3 , 4 , 5 , 6 , 1E8 // idea is to have the infinite volume guy here
} ;

#endif

/*
static const double L3[13] = {
  13824 , 32768 ,
  13824 , 32768 ,
  13824 , 32768 , 110592 ,
  110592 ,
  13824 , 32768 , 110592 , 262144 , 1E9 // idea is to have the infinite volume guy here
} ;
*/

static const double L3[13] = {
  0.008012086374634377  , 0.003314362727044057 ,
  0.008033356470775896  , 0.0033767579650790616 ,
  0.0033915369128406586 , 0.0010057346337476783 ,
  0.0010170795028159455 ,
  0.00813664 , 0.00343264 , 0.001017 , 0.0004291 , 0.0 , 0.0 // idea is to have the infinite volume guy here
} ;

double
ffvol1( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef V1
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars] ) ) ;
#elif (defined V2)
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars] )
			+ fparams[3]*exp( -sqrt(2)*MPIL[Npars] ) ) ;
#elif (defined THREEP1)
  return fparams[0]*X.X + fparams[1]*a2[Npars] ;
#elif (defined L3_EXTRAP) // 1/V volume dependence
    return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*L3[Npars] ) ;
#else
  return fparams[0] + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars] ) ;
#endif
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
    f[i] = ffvol1( X , p , i ) - DATA -> y[i] ;
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


#ifdef V1
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[i] ) ) ; 
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i] ) ;
#elif (defined V2)
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[2]*exp( -sqrt(2)*MPIL[i] ) ) ; 
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i] ) ;
    df[3][i] = fparams[0] * exp( -sqrt(2)*MPIL[i] ) ;
#elif (defined THREEP1)
    df[0][i] = X.X ;
    df[1][i] = a2[i] ;
#elif (defined L3_EXTRAP) // 1/V volume dependence
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[2]*L3[i] ) ; 
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * L3[i] ;
#else
    df[0][i] = 1 ;
    df[1][i] = X.X ;
    df[2][i] = exp( -MPIL[i] ) ;
#endif
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
  fparams[0] = 0.0 ;
  fparams[1] = 1 ;
  fparams[2] = 1 ;
  fparams[3] = 1 ;
  return ;
}
