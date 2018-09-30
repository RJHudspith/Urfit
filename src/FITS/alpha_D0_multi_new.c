/**
   @file alpha_D0.c
   @brief fit alpha_s from the D0 OPE


   Fit parameter map ::

   fparam[0] is the continuum coupling \alpha_s
   fparam[1] is the expected rotation-preserving lattice artifact -> fparam[1] * a^2(Q_1^2 - Q^2)/(t1-t)
   fparam[2] is the correction to alpha_s modelled as -> \alpha * ( 1 + fparam[3] * a^2 ) / Pi
   fparam[3] is the expected next-highest order correction -> fparam[2] * a^4(Q_1^4 - Q^4)/(t1-t)
   fparam[4] is a fit to the unknown parameter d_5 -> fparam[5] * [ \alpha * ( 1 + fparam[3] * a^2 ) / Pi ]^5

   It is a complicated fit.
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

const static double m[ 3*12 ] = { 2.0 , 2.0 , 2.0 , 2.25 , 2.25 , 2.25 , 2.5 , 2.5 , 2.5 , 2.75 , 2.75 , 2.75 ,
				  2.0 , 2.0 , 2.0 , 2.25 , 2.25 , 2.25 , 2.5 , 2.5 , 2.5 , 2.75 , 2.75 , 2.75 ,
				  2.0 , 2.0 , 2.0 , 2.25 , 2.25 , 2.25 , 2.5 , 2.5 , 2.5 , 2.75 , 2.75 , 2.75 } ;

const static double a2[ 3*12 ] = { 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 ,
				   0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 ,
				   0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 ,
				   0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 ,
				   0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 ,
				   0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 ,
				   0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 ,
				   0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 ,
				   0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 } ;

static double mu = 2.0 ;

static double Q1[ 3*12 ] = { 0 } ;

// beta parameters
static const double b0 = 2.25 , b1 = 4.0 , b2 = 10.059895833333334 , b3 = 47.228039589777325 ;

void
set_Q1_multi2( const double val , const size_t idx )
{
  Q1[ idx ] = val ;
  printf( "Q1 multi set :: %zu %f \n" , idx , Q1[idx] ) ;
  return ;
}

void
set_mu_multi2( const double munew )
{
  mu = munew ;
  return ;
}

const static double
lp1( const double t1 , const double t2 , const double d5 ) {
  return 1.0 ;
}

const static double
lp2( const double t1 , const double t2 , const double d5 ) {
  return 1.63982 - 1.125 * ( t1 + t2 ) ;
}

const static double
lp3( const double t1 , const double t2 , const double d5 ) {
  return 6.37108 - 5.6896 * ( t1 + t2 ) + 1.6875 * ( t1 * t1 + t1 * t2 + t2*t2 ) ;
}

const static double
lp4( const double t1 , const double t2 , const double d5 ) {
  return 49.0769 - 33.0914 * ( t1 + t2 ) + 15.8015 * ( t1 * t1 + t1 * t2 + t2*t2 ) \
    - 2.84766 * ( t1*t1*t1 + t1*t1*t2 + t2*t2*t1 + t2*t2*t2 ) ;
}

const static double
lp5( const double t1 , const double t2 , const double d5 ) {
  return d5 - 299.177 * ( t1 + t2 ) + 129.578 * ( t1*t1 + t1*t2 + t2*t2 ) \
    - 40.6161 * ( t1*t1*t1 + t1*t1*t2 + t2*t2*t1 + t2*t2*t2 )		\
    + 5.12578 * ( t1*t1*t1*t1 + t1*t1*t1*t2 + t1*t1*t2*t2 + t2*t2*t2*t1 + t2*t2*t2*t2 ) ;
}

const static double (*loop[5])( const double t1 , const double t2 , const double d5 ) = { lp1 , lp2 , lp3 , lp4 , lp5 } ;

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


double
falpha_D0_multi2( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double asq = 1.0 / ( fparams[4+Npars/12] * fparams[4+Npars/12] ) ; //a2[ Npars ] ;
  
  const double t1 = log( Q1[ Npars ] / ( m[ Npars ] * m[ Npars ] ) ) ;
  const double t2 = log( X.X / ( asq * m[ Npars ] * m[ Npars ] ) ) ;

  // Idea is to run alpha at mu to the other scale "m" and use that in the fit
  const double a_pi = rescale_alpha( fparams[0] / M_PI , mu , m[ Npars ] ) * ( 1 + fparams[2]*asq ) ;

  // momentum corrections
  double corrections = fparams[1] * ( X.X - Q1[ Npars ] * asq ) / ( t2 - t1 ) ;
  corrections += fparams[3] * ( X.X * X.X - Q1[ Npars] * Q1[ Npars]*asq*asq ) / ( t2 - t1 ) ;
  
  // delta from D=0 OPE
  register double PT = 0.0 ;
  size_t loops ;
  for( loops = LOOPS ; loops > 0 ; loops-- ) {
    PT = a_pi * ( loop[ loops-1 ](t1,t2,D5) + PT ) ;
  }  
  
  return PT + corrections ;
}

// D0_multi function evaluation
void
alpha_D0_multi2_f( double *f , const void *data , const double *fparams )
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
    f[i] = falpha_D0_multi2( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
alpha_D0_multi2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    size_t j ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( falpha_D0_multi2 , X ,
					   DATA -> map[i].bnd ,
					   fparams , j , DATA -> Npars ) ;
    }
 
  }
  
  return ;
}

void
alpha_D0_multi2_d2f( double **d2f , const void *data , const double *fparams )
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
alpha_D0_multi2_guesses( double *fparams ,
			const struct data_info Data ,
			const struct fit_info Fit )
{
  fparams[0] = 0.3 ;
  fparams[1] = 0.20 ;
  fparams[2] = -4.0 ;
  fparams[3] = -0.01 ;
  
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
