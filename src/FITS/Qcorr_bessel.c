/**
   @file Qcorr_Bessel.c
   @brief computes the topological correlator mass using

   y = A/x * K_1( mx )
   
   //y = m / ( 4\pi^2 * x) * e^{-mx} * \sqrt( \pi / (2mx) ) * ( 1 + 3/(8*mx))

   I am going to put the mass in fparams[0]
 */
#include "gens.h"

#include "gsl/gsl_sf_bessel.h"

#define FITAMP1

double
fQcorr_bessel( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef FITAMP1
  return fparams[0] / ( X.X ) * exp( -fparams[1] * X.X ) ;
#elif defined FITAMP2
  return fparams[0] / ( X.X ) * exp( -fparams[1] * X.X ) *
    sqrt( M_PI / ( 2. * fparams[1] * X.X ) )
    * ( 1 + 3./( 8 * fparams[1] * X.X ) )
    ;
#else
  return ( fparams[0] / ( 4.0 * M_PI * M_PI * X.X ) ) *
    exp( -fparams[0] * X.X ) *
    sqrt( M_PI / ( 2. * fparams[0] * X.X ) ) *
    ( 1 + 3./( 8 * fparams[0] * X.X ) ) ;
#endif
}

void
Qcorr_bessel_f( double *f , const void *data , const double *fparams )
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
    f[i] = fQcorr_bessel( X , p , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
Qcorr_bessel_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double h = 1E-8 ;
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      double p1[ DATA -> Npars ] , p2[ DATA -> Npars ] ;
      size_t k ;
      for( k = 0 ; k < DATA -> Npars ; k++ ) {
	p1[ k ] = p2[ k ] = fparams[ DATA -> map[ i ].p[ k ] ] ;
      }
      p1[j] += h ; p2[j] -= h ;
      
      df[ DATA -> map[i].p[ j ] ][i] =
	( fQcorr_bessel( X , p1 , DATA -> Npars ) -
	  fQcorr_bessel( X , p2 , DATA -> Npars ) ) / ( 2.*h ) ;
      
    }

  }
  return ;
}

void
Qcorr_bessel_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
Qcorr_bessel_guesses( double *fparams ,
	     const struct data_info Data ,
	     const struct fit_info Fit )
{
#if (defined FITAMP1) || (defined FITAMP2)
  fparams[0] = 0.02 ;
  fparams[1] = 0.3 ;
#else
  fparams[0] = 0.2 ;
#endif
}
