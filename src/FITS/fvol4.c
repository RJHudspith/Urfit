/**
   @file fvol4.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

#define PHIB

static const double MPIL[ 16 ] = {
  5.09, 4.00, 4.1056,
  5.83, 4.93, 3.88, 4.59,
  5.13, 4.32192 , 4.42, 3.802,
  6.44, 5.40, 4.43,
  5.07,
  1E5
} ;


// interesting \phi_B cancels a lot of the mpi-dependence in t0
#ifdef PHIB

static const double phiBcont = 120.084 ;

static const double PhiB[ 16 ] = {
  123.005, 123.508, 125.763,
  121.994, 122.646, 122.819, 123.249,
  123.420, 123.930, 123.210, 124.028,
  123.090, 122.900, 123.124,
  121.728,
  120.084
} ;

#else

// physb
static const double phiBcont = 5.2794 ;

static const double PhiB[ 15 ] = {
  5.287919646 , 5.273287378 , 5.25733375,
  5.287012435 , 5.279291420 , 5.266543405 , 5.255771915,
  5.30578488 , 5.28259226 ,
  5.313306005 , 5.306602012 , 5.302164417,
  5.28699184 ,
  5.2794
} ;

#endif

double
ffvol4( const struct x_desc X , const double *fparams , const size_t Npars )
{  
  return fparams[0] * ( 1 + fparams[1]*X.X
			+fparams[2]*exp( -MPIL[Npars] )
			+fparams[3]*fparams[5]*fparams[5]
			+fparams[4]*( PhiB[Npars] - phiBcont )) ;
}

void
fvol4_f( double *f , const void *data , const double *fparams )
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

    f[i] = ffvol4( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvol4_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  for( size_t i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    const double A   = fparams[ DATA -> map[ i ].p[ 0 ] ] ; 
    const double mpi = fparams[ DATA -> map[ i ].p[ 1 ] ] ;
    const double fv1 = fparams[ DATA -> map[ i ].p[ 2 ] ] ;
    const double B   = fparams[ DATA -> map[ i ].p[ 3 ] ] ; 
    const double phb = fparams[ DATA -> map[ i ].p[ 4 ] ] ;
    const double asq = fparams[ DATA -> map[ i ].p[ 5 ] ] * fparams[ DATA -> map[ i ].p[ 5 ] ] ;


    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*X.X + fv1*exp( -MPIL[i] ) + B*asq + phb*( PhiB[i] - phiBcont )) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * X.X ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MPIL[i] ) ;
    df[ DATA -> map[ i ].p[ 3 ] ][i] = A * asq ;
    df[ DATA -> map[ i ].p[ 4 ] ][i] = A*( PhiB[i] - phiBcont ) ;
    df[ DATA -> map[ i ].p[ 5 ] ][i] = 2* A * B * fparams[ DATA -> map[ i ].p[ 5 ] ] ;
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvol4_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvol4_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = -120.0 ;
  fparams[1] = 0,1 ;
  fparams[2] = 1 ;
  fparams[3] = 0 ;
  return ;
}
