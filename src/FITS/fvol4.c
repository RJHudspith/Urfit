/**
   @file fvol4.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

#define ORDERA
//#define PHIB

static const double MPIL[ 21 ] = {
  5.09, 4.00, 4.1056,
  5.83, 4.93, 3.88, 4.59,
  5.13, 4.32192 , 4.42, //3.802,
  6.44, 5.40, 4.43,
  5.07,
  1E5,1E5,1E5,1E5,1E5,1E5
} ;

static const double a[ 21 ] = {
  0.09929, 0.09929, 0.09929,
  0.08636, 0.08636, 0.08636, 0.08636,
  0.07634, 0.07634, 0.07634,
  0.06426, 0.06426, 0.06426,
  0.04981,
  0.09929 , 0.08636, 0.07634, 0.06426, 0.04981, 0.0,
} ;

// interesting \phi_B cancels a lot of the mpi-dependence in t0
#ifdef PHIB

static const double phiBcont = 120.084 ;

static const double PhiB[ 21 ] = {
  123.005, 123.508, 125.763,
  121.994, 122.646, 122.819, 123.249,
  123.420, 123.930, 123.210, //124.028,
  123.090, 122.900, 123.124,
  121.728,
  120.084,  120.084,  120.084,  120.084,  120.084,  120.084
} ;

#else

// physb
static const double phiBcont = 5.2794 ;

// redo this shit
static const double PhiB[ 21 ] = {
  5.285251105, 5.272395215, 5.25032958,
  5.287556265, 5.27939653, 5.266739915, 5.255630245,
  5.30536094, 5.296047185, 5.280573374,
  5.311340565, 5.306211995, 5.302109139,
  5.286052856,
  5.2794,  5.2794,  5.2794,  5.2794,  5.2794,  5.2794
} ;

#endif

double
ffvol4( const struct x_desc X , const double *fparams , const size_t Npars )
{  
  return fparams[0] * ( 1
			+ fparams[1]*X.X
			+fparams[2]*exp( -MPIL[Npars] )
			#ifdef ORDERA
			+fparams[3]*a[Npars] //fparams[5]
			#else
			+fparams[3]*a[Npars]*a[Npars] //fparams[5]*fparams[5]
			#endif
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

    #ifdef ORDERA
    const double asq = a[i] ; //fparams[ DATA -> map[ i ].p[ 5 ] ] ;
    #else
    const double asq = a[i]*a[i] ; //fparams[ DATA -> map[ i ].p[ 5 ] ] * fparams[ DATA -> map[ i ].p[ 5 ] ] ;
    #endif

    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*X.X + fv1*exp( -MPIL[i] ) + B*asq + phb*( PhiB[i] - phiBcont )) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * X.X ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MPIL[i] ) ;
    df[ DATA -> map[ i ].p[ 3 ] ][i] = A * asq ;
    df[ DATA -> map[ i ].p[ 4 ] ][i] = A*( PhiB[i] - phiBcont ) ;
    //df[ DATA -> map[ i ].p[ 5 ] ][i] = 2* A * B * fparams[ DATA -> map[ i ].p[ 5 ] ] ;
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
  fparams[1] = 0.1 ;
  fparams[2] = 1 ;
  fparams[3] = 0 ;
  return ;
}
