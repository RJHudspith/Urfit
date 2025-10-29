/**
   @file fvol4.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

static const double MPIL[ 15 ] = {
  3.5997872, 3.54, 3.4590432,
  3.442263, 3.5079968,
  3.704123, 3.4551776,
  3.476592,
  3.6751712,

  10000,19999,19999,19999,19999,19999
} ;

static const double Asq[ 15 ] = {
  0.36730945821854916, 0.36730945821854916, 0.36730945821854916,
  0.25767218944059367, 0.25767218944059367,
  0.1625910509885536, 0.1625910509885536,
  0.10078105316200554,
  0.06376900316294255,

  //
  //0,//0,0,
  0.36730945821854916, 0.25767218944059367, 0.1625910509885536, 0.10078105316200554, 0.06376900316294255, 0.0
} ;

static const double xcont = 0.07810406884989159 ;

double
ffvol6( const struct x_desc X , const double *fparams , const size_t Npars )
{  
  return fparams[0] * ( 1 + fparams[1]*(X.X-xcont)
			+fparams[2]*exp( -MPIL[Npars] ) //MPIL[Npars]
			+fparams[3]*Asq[Npars]
			+fparams[4]*Asq[Npars]*Asq[Npars]
			) ;
}

void
fvol6_f( double *f , const void *data , const double *fparams )
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

    f[i] = ffvol6( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvol6_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  for( size_t i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    const double A   = fparams[ DATA -> map[ i ].p[ 0 ] ] ; 
    const double mpi = fparams[ DATA -> map[ i ].p[ 1 ] ] ;
    const double fv1 = fparams[ DATA -> map[ i ].p[ 2 ] ] ;
    const double B   = fparams[ DATA -> map[ i ].p[ 3 ] ] ;
    const double C   = fparams[ DATA -> map[ i ].p[ 4 ] ] ;
    const double asq = Asq[i] ;
    const double a4  = Asq[i]*Asq[i] ;

    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*(X.X-xcont) + fv1*exp( -MPIL[i] )+ B*asq + C*a4 ) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * (X.X-xcont) ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MPIL[i] ) ; ///MPIL[i] ;
    df[ DATA -> map[ i ].p[ 3 ] ][i] = A * asq ;
    df[ DATA -> map[ i ].p[ 4 ] ][i] = A * a4 ;
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvol6_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvol6_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = -120.0 ;
  fparams[1] = 0,1 ;
  fparams[2] = 1 ;
  fparams[3] = 0 ;
  return ;
}
