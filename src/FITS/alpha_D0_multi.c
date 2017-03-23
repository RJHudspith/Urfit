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

#define LOOPS (5)

// alpha_s correction term == \alpha/pi * ( 1 + c_alpha a^2 )
#define ALPHA_CORR

// aQ^4 correction
#define FIT_A4

// aQ^6 correction
//#define FIT_A6

#define FIT_ACORR_A4

//#define FIT_D5

#ifndef FIT_D5
  #define D5 (400)
#endif

const static double m[ 3*12 ] = { 2.0 , 2.0 , 2.0 , 2.25 , 2.25 , 2.25 , 2.5 , 2.5 , 2.5 , 2.75 , 2.75 , 2.75 ,
				  2.0 , 2.0 , 2.0 , 2.25 , 2.25 , 2.25 , 2.5 , 2.5 , 2.5 , 2.75 , 2.75 , 2.75 ,
				  2.0 , 2.0 , 2.0 , 2.25 , 2.25 , 2.25 , 2.5 , 2.5 , 2.5 , 2.75 , 2.75 , 2.75 } ;

static double mu = 2.0 ;

static double Q1[ 3*12 ] = { 0 } ;
const static double a2[ 3*12 ] = { 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 ,
				   0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 ,
				   0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 , 0.10090915108763919 ,
				   0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 ,
				   0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 ,
				   0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 , 0.17605265301057807 ,
				   0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 ,
				   0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 ,
				   0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 , 0.31392137319354574 } ;

// beta parameters
const double b0 = 2.25 , b1 = 4.0 , b2 = 10.059895833333334 , b3 = 47.228039589777325 ;

void
set_Q1_multi( const double val , const size_t idx )
{
  Q1[ idx ] = val ;
  printf( "Q1 multi set :: %zu %f \n" , idx , Q1[idx] ) ;
  return ;
}

void
set_mu_multi( const double munew )
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
falpha_D0_multi( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double t1 = log( Q1[ Npars ] / ( m[ Npars ] * m[ Npars ] ) ) ;
  const double t2 = log( X.X / ( m[ Npars ] * m[ Npars ] ) ) ;

  // Idea is to run alpha at mu to the other scale "m" and use that in the fit
#ifdef ALPHA_CORR
  const double a_pi = rescale_alpha( fparams[0] / M_PI , mu , m[ Npars ] ) * ( 1 + a2[ Npars ] * fparams[2] ) ;
#else
  const double a_pi = rescale_alpha( fparams[0] / M_PI , mu , m[ Npars ] ) ;
#endif

  // momentum corrections
  double corrections = fparams[1] * a2[ Npars ] * ( Q1[ Npars] - X.X ) / ( t1 - t2 ) ;
#ifdef FIT_A4
  corrections += fparams[3] * a2[ Npars ] * a2[ Npars ] * ( Q1[ Npars] * Q1[ Npars] - X.X * X.X ) / ( t1 - t2 ) ;
#endif

#ifdef FIT_A6
  corrections += fparams[4] * a2[ Npars ] * a2[ Npars ] * ( Q1[ Npars] * Q1[ Npars] * Q1[ Npars ] - X.X * X.X * X.X ) / ( t1 - t2 ) ;
#endif
  
  // delta from D=0 OPE
  register double PT = 0.0 ;
  size_t loops ;
  for( loops = LOOPS ; loops > 0 ; loops-- ) {
    #ifdef FIT_D5
    PT = a_pi * ( loop[ loops-1 ](t1,t2,fparams[4]) + PT ) ;
    #else
    PT = a_pi * ( loop[ loops-1 ](t1,t2,D5) + PT ) ;
    #endif
  }  
  
  return PT + corrections ;
}

// D0_multi function evaluation
void
alpha_D0_multi_f( double *f , const void *data , const double *fparams )
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
    f[i] = falpha_D0_multi( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
alpha_D0_multi_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    // this is the data index
    const size_t mref = DATA -> map[i].bnd ;

    const double t1 = log( Q1[ mref ] / ( m[ mref ] * m[ mref ] ) ) ;
    const double t2 = log( DATA -> x[i] / ( m[ mref ] * m[ mref ] ) ) ;

    // cache of results for the fit
    const double asq = a2[ mref ] ;

#ifdef ALPHA_CORR
    const double acr = fparams[ DATA -> map[i].p[2] ] ;
    const double a_pi = rescale_alpha( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) * ( 1 + acr * asq ) ;
#else
    const double a_pi = rescale_alpha( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) ;
#endif

#ifdef FIT_D5
    const double fd5 = fparams[ DATA -> map[i].p[4] ] ;
#endif
    const double Qref = Q1[ mref ] ;

    // these two depend on the loop order of PT
    size_t loops ;
    register double dalpha = 0.0 , dacorr = 0.0 ;
    for( loops = LOOPS ; loops > 1 ; loops-- ) {
      #ifdef FIT_D5
      dalpha = a_pi * ( loops * loop[ loops-1 ](t1,t2,fd5) + dalpha ) ;
      #else
      dalpha = a_pi * ( loops * loop[ loops-1 ](t1,t2,D5) + dalpha ) ;
      #endif
    }
    dalpha += 1.0 ;
    
#ifdef ALPHA_CORR
    dacorr = rescale_alpha( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) * asq * dalpha ;
#endif
    dalpha = rescale_alpha_der( fparams[ DATA -> map[i].p[0] ] / M_PI , mu , m[ mref ] ) * dalpha ;

#ifdef ALPHA_CORR
    dalpha *= ( 1 + acr * asq ) ;
#endif

    // alpha_s and it's correction terms
    df[ DATA -> map[i].p[0] ][i] = dalpha ;

#ifdef ALPHA_CORR
    df[ DATA -> map[i].p[2] ][i] = dacorr ;
#endif
    
    // rotation-preserving corrections
    df[ DATA -> map[i].p[1] ][i] = asq *  ( Qref - DATA -> x[i] ) / ( t1 - t2 ) ;
#ifdef FIT_A4
    df[ DATA -> map[i].p[3] ][i] = asq * asq * ( Qref * Qref - DATA -> x[i] * DATA -> x[i] ) / ( t1 - t2 ) ;
#endif

#ifdef FIT_A6
    df[ DATA -> map[i].p[4] ][i] = asq * asq * asq * ( Qref * Qref * Qref - DATA -> x[i] * DATA -> x[i] * DATA -> x[i] ) / ( t1 - t2 ) ;
#endif
    
    // derivative of perturbative free parameter d_5
#ifdef FIT_D5
    df[ DATA -> map[i].p[4] ][i] = pow( a_pi , 5 ) ;
#endif
  }
  
  return ;
}

void
alpha_D0_multi_d2f( double **d2f , const void *data , const double *fparams )
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
alpha_D0_multi_guesses( double *fparams ,
		  const struct data_info Data ,
		  const struct fit_info Fit )
{
  fparams[0] = 0.3 ;
  fparams[1] = 0.17 ;
  
#ifdef ALPHA_CORR
  fparams[2] = -3.0 ;
#endif
  
#ifdef FIT_A4
  fparams[3] = -0.025 ;
#endif
  
#ifdef FIT_D5
  fparams[4] = 10000 ;
#endif

  return ;
}
