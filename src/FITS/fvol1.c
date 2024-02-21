/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

#define FVMPIL
//#define FVMKL
#define WU101
#define PHYSS
//#define NO_SU3
//#define HIGHER_MASS

#if (defined NO_SU3)
  #if (defined PHYSS)
    #if (defined WU101)
      #define DLEN (13)
    #else
      #define DLEN (12)
    #endif
      #elif (defined WU101)
    #define DLEN (11)
  #else
      #define DLEN (10)
  #endif
#else
#if (defined PHYSS)
  #if (defined WU101)
     #define DLEN (15)
  #else
     #define DLEN (14)
  #endif
#elif (defined WU101)
  #define DLEN (13)
#else
  #define DLEN (12)
#endif
#endif

// U103, H101, U102 , H102, H105, N101, C101
static const double MKL[ DLEN ] = {
#ifndef NO_SU3
  4.3269264 , 5.8321216 ,
#endif
  4.6145688 , 6.100112 ,
#ifdef WU101
  4.8274632 ,
#endif
  6.436992 , 9.6920208 ,
  9.8667504 ,
#ifdef PHYSS
  7.6133856 ,
  7.2080928 ,
#endif
  4.5 , 5.5 , 6.5 , 7.5 , 1E8  // idea is to have the infinite volume guy here
} ;

// U103, H101, U102 , H102, U101, H105, N101, C101, H107, H106
static const double MPIL[ DLEN ] = {
#ifndef NO_SU3
  4.3269264 , 5.8321216 ,
#endif
  3.7376688 , 4.9293248,
#ifdef WU101
  2.8801056 ,
#endif
  3.91064 , 5.8579536 ,
  4.6628376 ,
#ifdef PHYSS
  5.1151296 ,
  3.8789792 ,
#endif
  3 , 4 , 5 , 6 , 1E8  // idea is to have the infinite volume guy here
} ;

#ifdef HIGHER_MASS

#define SUB (0.07822077018599871)
#define SPOW (2)
const double sub[DLEN] =
{
  7.496458e-01 , 7.561123e-01 ,
  5.603588e-01 , 5.469004e-01 ,
  #ifdef WU101
  3.377914e-01 ,
  #endif
  3.450699e-01 , 3.444669e-01 ,
  2.194823e-01 ,
  #ifdef PHYSS
  0.5549746522603248 ,
  0.3305229618046193 ,
  #endif
  0 , 0 , 0 , 0 , 0
} ;
#endif

double
ffvol1( const struct x_desc X , const double *fparams , const size_t Npars )
{  
#ifdef FVMPIL
  #ifdef HIGHER_MASS
  //return fparams[0] * ( 1 + fparams[1]*X.X+ fparams[3]*X.X*X.X + fparams[2]*exp( -MPIL[Npars] ) ) ;
  return fparams[0] * ( 1 + fparams[1]*X.X
			+fparams[3]*(pow(sub[Npars] - SUB , SPOW ) )
			+fparams[2]*exp( -MPIL[Npars] ) ) ;
  #else
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars] ) ) ;
  #endif
#elif (defined FVMKL)
  #ifdef HIGHER_MASS
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[3]*X.X*X.X + fparams[2]*exp( -MKL[Npars] ) ) ;
  #else
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MKL[Npars] ) ) ;
  #endif
#else
  return fparams[0] * ( 1 + fparams[1]*X.X
			+fparams[2]*exp( -MPIL[Npars] )
			+fparams[3]*exp( -MKL[Npars] ) ) ;
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

    const size_t mi = i%(DLEN-5) ;
    f[i] = ffvol1( X , p , mi ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvol1_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  for( size_t i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    const double A   = fparams[ DATA -> map[ i ].p[ 0 ] ] ; 
    const double mpi = fparams[ DATA -> map[ i ].p[ 1 ] ] ;
    const double fv1 = fparams[ DATA -> map[ i ].p[ 2 ] ] ; 

    const size_t mi = i%(DLEN-5) ;
    
#if (defined FVMPIL)
    #ifdef HIGHER_MASS
    const double mpisq = fparams[ DATA -> map[ i ].p[ 3 ] ] ;
    const double mterm = (pow(sub[ mi ] - SUB , SPOW ) ) ;
    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*X.X
					 + mpisq*mterm
					 + fv1*exp( -MPIL[mi] ) ) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * X.X ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MPIL[mi] ) ;
    df[ DATA -> map[ i ].p[ 3 ] ][i] = A * mterm ;
    #else
    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*X.X + fv1*exp( -MPIL[mi] ) ) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * X.X ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MPIL[mi] ) ;
    #endif
#elif (defined FVMKL)
    #ifdef HIGHER_MASS
    const double mpisq = fparams[ DATA -> map[ i ].p[ 3 ] ] ;
    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*X.X + mpisq*X.X*X.X + fv1*exp( -MPIL[mi] ) ) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * X.X ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MKL[mi] ) ;
    df[ DATA -> map[ i ].p[ 3 ] ][i] = A * X.X*X.X ;
    #else
    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*X.X + fv1*exp( -MKL[mi] ) ) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * X.X ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MKL[mi] ) ;
    #endif
#else
    const double fv2 = fparams[ DATA -> map[ i ].p[ 3 ] ] ; 
    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*X.X
		 + fv1*exp( -MPIL[mi] )
		 + fv2*exp( -MKL[mi] ) ) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * X.X ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MPIL[mi] ) ;
    df[ DATA -> map[ i ].p[ 3 ] ][i] = A * exp( -MKL[mi] ) ;
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
  //fparams[3] = 1 ;
  return ;
}
