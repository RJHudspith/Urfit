/**
   @file alpha_D0.c
   @brief fit alpha_s from the D0 OPE
 */
#include "gens.h"

#include "Nder.h"
#include "cruel_runnings.h"

#define LOOPS (4)

// alpha_s correction term == \alpha/pi * ( 1 + c_alpha a^2 )
#define ALPHA_CORR

// aQ^4 correction
#define FIT_AQ4

#ifndef FIT_D5
  #define D5 (400)
#endif

static double mu = 2.0 ;

static const double a2[ 9 ] = { 0.04949440885942718 , 0.04949440885942718 , 0.04949440885942718 ,
				0.07684254741528086 , 0.07684254741528086 , 0.07684254741528086 ,
				0.16601189178800163 , 0.16601189178800163 , 0.16601189178800163 } ;

const static double m[ 9 ] = { 2.9 , 3 , 3.1 ,
			       2.9 , 3 , 3.1 ,
			       2.9 , 3 , 3.1 } ;

// beta parameters
static const double b0 = 2.25 , b1 = 4.0 , b2 = 10.059895833333334 , b3 = 47.228039589777325 ;

const static double
lp1( const double t , const double d5 ) {
  return 1.0 ;
}

const static double
lp2( const double t , const double d5 ) {
  return 1.63982 - 2.25 * t ;
}

const static double
lp3( const double t , const double d5 ) {
  return 6.37108 - 11.379195 * t + 5.0625 * t ;
}

const static double
lp4( const double t , const double d5 ) {
  return 49.0769 -66.182813*t + 47.404784*t*t -11.390625*t*t*t ;
}

const static double
lp5( const double t , const double d5 ) {
  return d5 -598.354375*t + 388.732597*t*t - 162.464353*t*t*t + 25.628906*t*t*t*t ;
}

const static double (*loop[5])( const double t , const double d5 ) = { lp1 , lp2 , lp3 , lp4 , lp5 } ;

// writes alpha(m) in terms of alpha(\mu) == a_pim
static double
rescale_alpha( const double a_pim , const double m , const double mup )
{
  const double t = log( mup*mup/(m*m) ) ;
#if LOOPS==1
  return a_pim * ( 1 ) ;
#elif LOOPS==2
  return a_pim * ( 1 + a_pim * ( -b0 * t ) ) ;
#elif LOOPS==3
  return a_pim * ( 1 + a_pim * ( -b0 * t + a_pim * ( ( b0*b0*t*t - b1 * t ) ) ) ) ;
#elif LOOPS==4
  return a_pim * ( 1 + a_pim * ( -b0 * t + a_pim * ( ( b0*b0*t*t - b1 * t ) + a_pim * ( ( -b0*b0*b0*t*t*t + 5/2.*b0*b1*t*t - b2 *t ) ) ) ) ) ;
#elif LOOPS==5
  return a_pim * ( 1 + a_pim * ( -b0 * t + a_pim * ( ( b0*b0*t*t - b1 * t ) + a_pim * ( ( -b0*b0*b0*t*t*t + 5/2.*b0*b1*t*t - b2 *t ) + a_pim * ( b0*b0*b0*b0*t*t*t*t - 13/3.*b0*b0*b1*t*t*t + (3*b1*b1+6*b0*b2)/2.*t*t - b3*t ) ) ) ) ) ;
#endif
}

// derivative wrt to alpha(\mu) of \alpha(m)
static double
rescale_alpha_der( const double a_pim , const double m , const double mup )
{
  const double t = log( mup*mup/(m*m) ) ;
#if LOOPS==1
  return ( 1 ) / M_PI ;
#elif LOOPS==2
  return ( 1 + a_pim * ( -2 * b0 * t ) ) / M_PI ;
#elif LOOPS==3
  return ( 1 + a_pim * ( -2 * b0 * t + a_pim * ( 3 * ( b0*b0*t*t - b1 * t ) ) ) ) / M_PI ;
#elif LOOPS==4
  return ( 1 + a_pim * ( -2 * b0 * t + a_pim * ( 3 * ( b0*b0*t*t - b1 * t ) + a_pim * ( 4 * ( -b0*b0*b0*t*t*t + 5/2.*b0*b1*t*t - b2 *t ) ) ) ) ) / M_PI ;
#elif LOOPS==5
  return ( 1 + a_pim * ( -2 * b0 * t + a_pim * ( 3 * ( b0*b0*t*t - b1 * t ) + a_pim * ( 4 * ( -b0*b0*b0*t*t*t + 5/2.*b0*b1*t*t - b2 *t ) + a_pim * 5 * ( b0*b0*b0*b0*t*t*t*t - 13/3.*b0*b0*b1*t*t*t + (3*b1*b1+6*b0*b2)/2.*t*t - b3*t ) ) ) ) ) / M_PI ;
#endif
}

void
set_mu_multi_adler( const double munew )
{
  mu = munew ;
  return ;
}

double
fadleralpha_D0_multi( const struct x_desc X , const double *fparams , const size_t Npars )
{
  double a2inv = fparams[4+Npars/3] * fparams[4+Npars/3] ;

  //printf( "A2INV %zu %e \n" , Npars , sqrt( a2inv) ) ;
  
  const double t = log( a2inv * X.X / ( m[ Npars ] * m[ Npars ] ) ) ;
  
  // Idea is to run alpha at mu to the other scale "m" and use that in the fit
  const double a_pi = rescale_alpha( fparams[0] / M_PI , mu , m[ Npars ] ) ;
  
  // momentum corrections
  double corrections =
    fparams[3] / a2inv +
    fparams[1] * ( X.X ) +
    fparams[2] * ( X.X * X.X ) +
    fparams[10] * ( X.X * X.X * X.X ) ;
 
  // delta from D=0 OPE
  register double PT = 0.0 ;
  size_t loops ;
  for( loops = LOOPS ; loops > 0 ; loops-- ) {
    PT = a_pi * ( loop[ loops-1 ](t,D5) + PT ) ;
  }

  //const double ZV = ( 1. - fparams[7] / a2inv ) ;
  const double ZV = fparams[7+Npars/3] ;
  return ZV * ZV * ( PT + corrections ) ;
}

// D0_multi function evaluation
void
adleralpha_D0_multi_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> Npars ] ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      p[j] = fparams[ DATA -> map[i].p[j] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fadleralpha_D0_multi( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
adleralpha_D0_multi_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    size_t j ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( fadleralpha_D0_multi , X ,
					   DATA -> map[i].bnd ,
					   fparams , j , DATA -> Npars ) ;
    }

#if 0

    // this is the data index
    const size_t mref = DATA -> map[i].bnd ;

    const double t = log( DATA -> x[i] / ( m[ mref ] * m[ mref ] ) ) ;

    // cache of results for the fit
    const double asq = a2[ mref ] ;

    const double a_pi = rescale_alpha( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) ;

    // these two depend on the loop order of PT
    size_t loops ;
    register double dalpha = 0.0 ;
    for( loops = LOOPS ; loops > 1 ; loops-- ) {
      dalpha = a_pi * ( loops * loop[ loops-1 ](t,D5) + dalpha ) ;
    }
    dalpha += 1.0 ;
    
    // alpha_s and its correction terms
    df[ DATA -> map[i].p[0] ][i] = rescale_alpha_der( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) * dalpha ;
    
    // rotation-preserving corrections
    df[ DATA -> map[i].p[1] ][i] = asq *  ( DATA -> x[i] ) ;

    df[ DATA -> map[i].p[2] ][i] = asq ;
    
    df[ DATA -> map[i].p[3] ][i] = asq * asq * ( DATA -> x[i] * DATA -> x[i] ) ;

    // set dalpha to zero
    dalpha = 0.0 ;
    for( loops = LOOPS ; loops > 0 ; loops-- ) {
      dalpha = a_pi * ( loops * loop[ loops-1 ](t,D5) + dalpha ) ;
    }
    df[ 4 + mref/3 ][i] = 2 * fparams[4 + mref/3 ] * dalpha ;

    /*
    // this is the data index
    const size_t mref = DATA -> map[i].bnd ;

    const double t = log( DATA -> x[i] / ( m[ mref ] * m[ mref ] ) ) ;

    // cache of results for the fit
    const double asq = a2[ mref ] ;

    const double acr = fparams[ DATA -> map[i].p[2] ] ;
    const double a_pi = rescale_alpha( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) * ( 1 + acr * asq ) ;

    // these two depend on the loop order of PT
    size_t loops ;
    register double dalpha = 0.0 , dacorr = 0.0 ;
    for( loops = LOOPS ; loops > 1 ; loops-- ) {
      dalpha = a_pi * ( loops * loop[ loops-1 ](t,D5) + dalpha ) ;
    }
    dalpha += 1.0 ;
    
    dacorr = rescale_alpha( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) * asq * dalpha ;

    dalpha *= ( 1 + acr * asq ) ;

    // alpha_s and it's correction terms
    df[ DATA -> map[i].p[0] ][i] = dalpha ;

    df[ DATA -> map[i].p[2] ][i] = dacorr ;
    
    // rotation-preserving corrections
    df[ DATA -> map[i].p[1] ][i] = asq *  ( DATA -> x[i] ) ;
    df[ DATA -> map[i].p[3] ][i] = asq * asq * ( DATA -> x[i] * DATA -> x[i] ) ;
    */
#endif
  }
  
  return ;
}

void
adleralpha_D0_multi_d2f( double **d2f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    for( j = 0 ; j < DATA -> Npars * DATA -> Npars ; j++ ) {
      d2f[j][i] = 0.0 ;
    }
  }
  return ;
}

void
adleralpha_D0_multi_guesses( double *fparams ,
			const struct data_info Data ,
			const struct fit_info Fit )
{
  fparams[0] = 0.3 ;
  fparams[1] = -0.20 ;
  fparams[2] = -2.5 ;
  fparams[3] = -0.01 ;
  fparams[4] = 1 ;

  // if we have some priors we set them
  size_t i ;
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    if( Fit.Prior[i].Initialised == true ) {
      fparams[i] = Fit.Prior[i].Val ;
    }
  }

  // print guesses
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    fprintf( stdout , "[ALPHA_D0_MULTI] Guesses :: %zu %f \n" ,
	     i , fparams[i] ) ;
  }

  return ;
}
