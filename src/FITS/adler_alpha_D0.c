/**
   @file alpha_D0.c
   @brief fit alpha_s from the D0 adler function
 */
#include "gens.h"

#include "Nder.h"

#define LOOPS (4)

#ifndef FIT_D5
  #define D5 (400)
#endif

static double mu = 2.0 ;
static double Q1 = 4.0 ;

static const double a2[ 3 ] = { 0.04949440885942718 ,
				0.07684254741528086 ,
				0.16601189178800163 } ;

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

void
set_mu_adleralpha( const double munew )
{
  mu = munew ;
  printf( "Mu set %f \n" , mu ) ;
  return ;
}

void
set_Q1_adleralpha( const double Q1new )
{
  Q1 = Q1new ;
  printf( "Q1 set %f \n" , Q1 ) ;
  return ;
}

double
fadleralpha_D0( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double t = log( X.X / ( mu * mu ) ) ;

  // program in an a_pi correction?
  double a_pi = fparams[0] / M_PI ;
  const double corrections = fparams[1] * a2[ Npars ] * ( X.X ) + 
    fparams[2] * a2[ Npars ] * a2[ Npars ] * ( X.X * X.X ) ;
  
  // delta from D=0 OPE  
  register double PT = 0.0 ;
  size_t loops ;
  for( loops = LOOPS ; loops > 0 ; loops-- ) {
    PT = a_pi * ( loop[ loops-1 ](t,D5) + PT ) ;
  }
  
  return ( PT + corrections + a2[ Npars ] * fparams[3] ) ;
}

// D0 function evaluation
void
adleralpha_D0_f( double *f , const void *data , const double *fparams )
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
    f[i] = fadleralpha_D0( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
adleralpha_D0_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    size_t j ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( fadleralpha_D0 , X ,
					   DATA -> map[i].bnd ,
					   fparams , j , DATA -> Npars ) ;
    }


#if 0
    const double t = log( DATA -> x[i] / ( mu * mu ) ) ;

    // cache of results for the fit
    const double asq = a2[ DATA -> map[i].bnd ] ;
    const double acr = fparams[ DATA -> map[i].p[3] ] ;
    const double a_pi = fparams[ DATA -> map[i].p[0] ] *  ( 1 + acr * asq ) / M_PI ;
    //const double fd5 = fparams[ DATA -> map[i].p[4] ] ;

    // these two depend on the loop order of PT
    size_t loops ;
    register double dalpha = 0.0 , dacorr = 0.0 ;
    for( loops = LOOPS ; loops > 1 ; loops-- ) {
      #ifdef FIT_D5
      dalpha = a_pi * ( loops * loop[ loops-1 ](t,fparams[4]) + dalpha ) ;
      #else
      dalpha = a_pi * ( loops * loop[ loops-1 ](t,D5) + dalpha ) ;
      #endif
    }
    #ifdef FIT_D5
    dacorr = fparams[ DATA -> map[i].p[0] ] * asq * ( loop[0](t,fparams[4]) + dalpha ) / M_PI ;
    dalpha = ( 1 + acr * asq ) * ( loop[ 0 ](t,fparams[4]) + dalpha ) / M_PI ;
    #else
    dacorr = fparams[ DATA -> map[i].p[0] ] * asq * ( loop[0](t,D5) + dalpha ) / M_PI ;
    dalpha = ( 1 + acr * asq ) * ( loop[ 0 ](t,D5) + dalpha ) / M_PI ;
    #endif

    // alpha_s and it's correction terms
    df[ DATA -> map[i].p[0] ][i] = dalpha ;
    df[ DATA -> map[i].p[3] ][i] = dacorr ;
    
    // rotation-preserving corrections
    df[ DATA -> map[i].p[1] ][i] = asq *  ( DATA -> x[i] ) ;
    df[ DATA -> map[i].p[2] ][i] = asq * asq * ( DATA -> x[i] * DATA -> x[i] ) ;

    // derivative of perturbative free parameter d_5
#ifdef FIT_D5
    df[ DATA -> map[i].p[4] ][i] = pow( a_pi , 5 ) ;
#endif
#endif
  }
  return ;
}

void
adleralpha_D0_d2f( double **d2f , const void *data , const double *fparams )
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
adleralpha_D0_guesses( double *fparams ,
		       const struct data_info Data ,
		       const struct fit_info Fit )
{
  fparams[0] = 0.3 ;
  fparams[1] = 0.17 ;
  fparams[2] = -0.019 ;
  fparams[3] = -3.0 ;
  return ;
}
