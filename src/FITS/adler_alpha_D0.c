/**
   @file alpha_D0.c
   @brief fit alpha_s from the D0 adler function
 */
#include "gens.h"

#include "Nder.h"

#define LOOPS (5)

#ifndef FIT_D5
  #define D5 (400)
#endif

static double mu = 2.0 ;

static const double a2[ 3 ] = { 0.04949440885942718 ,
				0.07684254741528086 ,
				0.16601189178800163 } ;

static const double
lp1( const double t , const double d5 ) {
  return 1.0 ;
}

static const double
lp2( const double t , const double d5 ) {
  return 1.63982 - 2.25 * t ;
}

static const double
lp3( const double t , const double d5 ) {
  return 6.37108 - 11.379195 * t + 5.0625 * t ;
}

static const double
lp4( const double t , const double d5 ) {
  return 49.0769 -66.182813*t + 47.404784*t*t -11.390625*t*t*t ;
}

static const double
lp5( const double t , const double d5 ) {
  return d5 -598.354375*t + 388.732597*t*t - 162.464353*t*t*t + 25.628906*t*t*t*t ;
}

static const double (*loop[5])( const double t , const double d5 ) = { lp1 , lp2 , lp3 , lp4 , lp5 } ;

void
set_mu_adleralpha( const double munew )
{
  mu = munew ;
  printf( "Mu set %f \n" , mu ) ;
  return ;
}

double
fadleralpha_D0( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double t = log( X.X / ( mu * mu ) ) ;

  const double a_pi = fparams[0] * ( 1. ) / M_PI ;
#ifdef BEST
  // program in an a_pi correction?
  const double corrections =
    fparams[1] * ( a2[Npars] * sqrt(a2[Npars]) *  X.X ) +
    fparams[3] * a2[ Npars ] * a2[Npars] * X.X * ( X.X )
    + fparams[2] * sqrt(a2[ Npars ]) * log( sqrt(a2[ Npars ]) ) 
    ;
#else
  // program in an a_pi correction?
  const double corrections =
    fparams[1] * ( a2[Npars] * sqrt(a2[Npars]) *  X.X ) +
    fparams[3] * a2[ Npars ] * a2[Npars] * X.X * ( X.X ) +
    fparams[2] * sqrt(a2[ Npars ]) * log( sqrt(a2[ Npars ]) * mu )
    //fparams[2] * a2[Npars]* a2[ Npars ] //log( sqrt(a2[ Npars ]) ) 
    ;
#endif
  
  // delta from D=0 OPE  
  register double PT = 0.0 ;
  size_t loops ;
  for( loops = LOOPS ; loops > 0 ; loops-- ) {
    PT = a_pi * ( loop[ loops-1 ](t,D5) + PT ) ;
  }
  
  return ( PT + corrections ) ;
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
