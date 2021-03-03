/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

#define POLE
//#define LOGSQ
//#define LOG
//#define MPI2LOGMPI2

// A653 , A654
// U103 , H101 , U102 , U101 , H105 , C101
// B450 , D450
// H200 , N202 , N200 , D200
// N300
static const double MPIL[16] = { 5.312 , 4.029 ,
				 4.354 , 5.818 , 3.743 , 3.917 , 4.642 ,
				 5.146 , 5.377 ,
				 4.358 , 6.408 , 4.416 , 4.154 ,
				 5.109 ,
				 1E8 ,
				 1E8 } ;

static const double a2[16] = { 0.25319 , 0.25319 ,
			       0.19536 , 0.19536 , 0.19536 , 0.19536 , 0.19536 ,
			       0.14967 , 0.14967 ,
			       0.10605 , 0.10605 , 0.10605 , 0.10605 ,
			       0.06372 ,
			       0.19536 ,
			       0.0 } ;

/*
static const double MPIL[11] = { 5.312 , 4.029 ,
				 4.354 , 5.818 , 3.917 , 4.642 ,
				 5.146 ,
				 4.358 , 6.408 ,
				 5.109 ,
				 1E8 } ;

static const double a2[11] = { 0.25319 , 0.25319 ,
			       0.19536 , 0.19536 , 0.19536 , 0.19536 ,
			       0.14967 ,
			       0.10605 , 0.10605 ,
			       0.06372 ,
			       0.0 } ;
*/

double
ffvol3( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef POLE
  return ( fparams[0]/X.X + fparams[1] + fparams[2]*X.X ) *
    ( 1 + fparams[3]*exp( -MPIL[Npars]/2 ) + fparams[4]*a2[Npars] ) ;
#elif (defined LOGSQ)
  return ( fparams[0]*log(X.X)*log(X.X) + fparams[1] + fparams[2]*X.X ) *
    ( 1 + fparams[3]*exp( -MPIL[Npars]/2 ) + fparams[4]*a2[Npars] ) ;
#elif (defined LOG)
  return ( fparams[0]*log(X.X) + fparams[1] + fparams[2]*X.X ) *
    ( 1 + fparams[3]*exp( -MPIL[Npars]/2 ) + fparams[4]*a2[Npars] ) ;
#elif (defined  MPI2LOGMPI2)
  return ( fparams[0]*X.X*log(X.X) + fparams[1] + fparams[2]*X.X ) *
    ( 1 + fparams[3]*exp( -MPIL[Npars]/2 ) + fparams[4]*a2[Npars] ) ;
#else
  return fparams[0]/X.X + fparams[1] + fparams[2]*X.X
    + fparams[3]*exp( -MPIL[Npars]/2 ) + fparams[4]*a2[Npars] ;
#endif
}

void
fvol3_f( double *f , const void *data , const double *fparams )
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
    f[i] = ffvol3( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvol3_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

#if (defined POLE)

    const double brac1 = ( fparams[0]/X.X + fparams[1] + fparams[2]*X.X ) ;
    const double brac2 = ( 1 + fparams[3]*exp( -MPIL[i]/2 ) + fparams[4]*a2[i] ) ;

    df[0][i] = brac2/X.X ;
    df[1][i] = brac2 ;
    df[2][i] = brac2*X.X ;

    df[3][i] = brac1*exp( -MPIL[i]/2 ) ;
    df[4][i] = brac1*a2[i] ;
    
#elif (defined LOGSQ)

    const double brac1 = ( fparams[0]*log(X.X)*log(X.X) + fparams[1] + fparams[2]*X.X ) ;
    const double brac2 = ( 1 + fparams[3]*exp( -MPIL[i]/2 ) + fparams[4]*a2[i] ) ;

    df[0][i] = brac2*log(X.X)*log(X.X);
    df[1][i] = brac2 ;
    df[2][i] = brac2*X.X ;

    df[3][i] = brac1*exp( -MPIL[i]/2 ) ;
    df[4][i] = brac1*a2[i] ;

#elif (defined LOG)

    const double brac1 = ( fparams[0]*log(X.X) + fparams[1] + fparams[2]*X.X ) ;
    const double brac2 = ( 1 + fparams[3]*exp( -MPIL[i]/2 ) + fparams[4]*a2[i] ) ;

    df[0][i] = brac2*log(X.X) ;
    df[1][i] = brac2 ;
    df[2][i] = brac2*X.X ;

    df[3][i] = brac1*exp( -MPIL[i]/2 ) ;
    df[4][i] = brac1*a2[i] ;

#elif (defined MPI2LOGMPI2)

    const double brac1 = ( fparams[0]*X.X*log(X.X) + fparams[1] + fparams[2]*X.X ) ;
    const double brac2 = ( 1 + fparams[3]*exp( -MPIL[i]/2 ) + fparams[4]*a2[i] ) ;

    df[0][i] = brac2*X.X*log(X.X) ;
    df[1][i] = brac2 ;
    df[2][i] = brac2*X.X ;

    df[3][i] = brac1*exp( -MPIL[i]/2 ) ;
    df[4][i] = brac1*a2[i] ;
    
    
#elif (defined V3)
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[4]/X.X + fparams[2]*exp( -MPIL[i]/2 ) + fparams[3]*a2[i] ) ;
    df[1][i] = fparams[0]*X.X ;
    df[2][i] = fparams[0]*exp( -MPIL[i]/2 ) ;
    df[3][i] = fparams[0]*a2[i] ;
    df[4][i] = fparams[0]/X.X ;
#else
    df[0][i] = 1/X.X ;
    df[1][i] = 1 ;
    df[2][i] = X.X ;

    df[3][i] = exp( -MPIL[i]/2 ) ;
    df[4][i] = a2[i] ;
#endif
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvol3_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvol3_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 10 ;
  fparams[1] = 80 ;
  fparams[2] = 10 ;
  fparams[3] = -0.1 ;
  fparams[4] = -0.1 ;
  
  return ;
}
