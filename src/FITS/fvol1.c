/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

// A653 , U103 , H101 , B450 , H200 , N202 , N300

static const double MPIL[8] = {  5.312 ,
				 4.354 , 5.818 ,
				 5.146 ,
				 4.358 , 6.408 ,
				 5.109 ,
				 1E8 // idea is to have the infinite volume guy here
} ;

static const double a2[9] = {
  0.25319 ,
  0.19536 , 0.19536 , 0.19536 ,
  0.14967 ,
  0.10605 , 0.10605 , 0.10605 ,
  0.0 } ;

#define THREEP1

double
ffvol1( const struct x_desc X , const double *fparams , const size_t Npars )
{
  //return fparams[0] + fparams[1]*X.X + fparams[2]*exp(-pow(MPIL[Npars],MPOW))/sqrt(MPIL[Npars]) ;
#ifdef V1
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars]/2 ) ) ;
#elif (defined V2)
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars]/2 )
			+ fparams[3]*exp( -sqrt(2)*MPIL[Npars]/2 ) ) ;
#elif (defined THREEP1)
  return fparams[0]*X.X + fparams[1]*a2[Npars] ;
#else
  return fparams[0] + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars]/2 ) ;
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
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[i]/2 ) ) ; 
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i]/2 ) ;
#elif (defined V2)
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[2]*exp( -sqrt(2)*MPIL[i]/2 ) ) ; 
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i]/2 ) ;
    df[3][i] = fparams[0] * exp( -sqrt(2)*MPIL[i]/2 ) ;
#elif (defined THREEP1)
    df[0][i] = X.X ;
    df[1][i] = a2[i] ;
    #else
    df[0][i] = 1 ; //( 1 + fparams[1] * X.X + fparams[2] *  exp( -MPIL[i]/2. ) ) ;
    df[1][i] = X.X ;
    df[2][i] = exp( -MPIL[i]/2 ) ;
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
