/**
   @file fvol4.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

#define ASQ

//#define WITH_COARSE

static const double MPIL[ 7 ] = {
#ifdef WITH_COARSE
  3.925310162857608,
#endif
  4.553798776021403,
#ifdef WITH_LOW_MPIL
  3.2893626608991835,
#endif
  4.37582088369264, 5.485644183071164,
  4.507345118239786, 4.791676434717644
} ;

static const double Aval[ 7 ] = {
#ifdef WITH_COARSE
  0.151,
#endif
  0.1207,
#ifdef WITH_LOW_MPIL
  0.1202,
#endif
  0.1184, 0.1189,
  0.0888, 0.0872
} ;

static const double xcont = 0.135*0.135 ;

double
ffvol6( const struct x_desc X , const double *fparams , const size_t Npars )
{  
  return fparams[0] * ( 1 + fparams[1]*(X.X-xcont)
			+fparams[2]*exp( -MPIL[Npars] ) //MPIL[Npars]
			#ifdef ASQ
			+fparams[3]*Aval[Npars]*Aval[Npars]
			#else
			+fparams[3]*Aval[Npars]
			#endif
			+fparams[4]*Aval[Npars]*Aval[Npars]*Aval[Npars]*Aval[Npars]
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
    #ifdef ASQ
    const double asq = Aval[i]*Aval[i] ;
    #else
    const double asq = Aval[i] ;
    #endif
    const double a4 = Aval[i]*Aval[i]*Aval[i]*Aval[i] ;

    df[ DATA -> map[ i ].p[ 0 ] ][i] = ( 1 + mpi*(X.X-xcont) + fv1*exp( -MPIL[i] )+ B*asq + C*a4 ) ; 
    df[ DATA -> map[ i ].p[ 1 ] ][i] = A * (X.X-xcont) ;
    df[ DATA -> map[ i ].p[ 2 ] ][i] = A * exp( -MPIL[i] )/MPIL[i] ;
    df[ DATA -> map[ i ].p[ 3 ] ][i] = A * asq ;
    df[ DATA -> map[ i ].p[ 4 ] ][i] = A * a4 ;
    /*
    #ifdef ASQ
    df[ DATA -> map[ i ].p[ 4 ] ][i] = A * B * 2* fparams[ DATA -> map[ i ].p[ 4 ] ] ;
    #else
    df[ DATA -> map[ i ].p[ 4 ] ][i] = A * B ;
    #endif
    */
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
