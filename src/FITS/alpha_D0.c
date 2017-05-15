/**
   @file alpha_D0.c
   @brief fit alpha_s from the D0 OPE


   Fit parameter map ::

   fparam[0] is the continuum coupling \alpha_s
   fparam[1] is the expected rotation-preserving lattice artifact -> fparam[1] * a^2(Q_1^2 - Q^2)/(t1-t)
   fparam[2] is the expected next-highest order correction -> fparam[2] * a^4(Q_1^4 - Q^4)/(t1-t)
   fparam[3] is the correction to alpha_s modelled as -> \alpha * ( 1 + fparam[3] * a^2 ) / Pi
   fparam[4] is a fit to the unknown parameter d_5 -> fparam[5] * [ \alpha * ( 1 + fparam[3] * a^2 ) / Pi ]^5

   It is a complicated fit.
 */
#include "gens.h"

#define LOOPS (5)

#ifndef FIT_D5
  #define D5 (400)
#endif

static double mu = 2.0 ;
static double Q1[ 3 ] ;
const static double a2[3] = { 0.10090915108763919 ,
			      0.17605265301057807 ,
			      0.31392137319354574 } ;

void
set_mu( const double munew )
{
  mu = munew ;
  printf( "Mu set %f \n" , mu ) ;
  return ;
}

void
set_Q1( const double val , const size_t idx )
{
  Q1[ idx ] = val ;
  printf( "Q1 set :: %zu %f \n" , idx , Q1[idx] ) ;
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

double
falpha_D0( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double t1 = log( Q1[ Npars ] / ( mu * mu ) ) ;
  const double t2 = log( X.X / ( mu * mu ) ) ;

  // program in an a_pi correction?
  double a_pi = fparams[0] / M_PI * ( 1 + a2[ Npars ] * fparams[3] ) ;
  double corrections = 0 ;

  // momentum corrections
  corrections += fparams[1] * a2[ Npars ] * ( Q1[ Npars] - X.X ) / ( t1 - t2 ) ;
  corrections += fparams[2] * a2[ Npars ] * a2[ Npars ] * ( Q1[ Npars] * Q1[ Npars] - X.X * X.X ) / ( t1 - t2 ) ;
  
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

// D0 function evaluation
void
alpha_D0_f( double *f , const void *data , const double *fparams )
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
    f[i] = falpha_D0( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
alpha_D0_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double t1 = log( Q1[ DATA -> map[i].bnd ] / ( mu * mu ) ) ;
    const double t2 = log( DATA -> x[i] / ( mu * mu ) ) ;

    // cache of results for the fit
    const double asq = a2[ DATA -> map[i].bnd ] ;
    const double acr = fparams[ DATA -> map[i].p[3] ] ;
    const double a_pi = fparams[ DATA -> map[i].p[0] ] *  ( 1 + acr * asq ) / M_PI ;
    //const double fd5 = fparams[ DATA -> map[i].p[4] ] ;
    const double Qref = Q1[ DATA -> map[i].bnd ] ;

    // these two depend on the loop order of PT
    size_t loops ;
    register double dalpha = 0.0 , dacorr = 0.0 ;
    for( loops = LOOPS ; loops > 1 ; loops-- ) {
      #ifdef FIT_D5
      dalpha = a_pi * ( loops * loop[ loops-1 ](t1,t2,fparams[4]) + dalpha ) ;
      #else
      dalpha = a_pi * ( loops * loop[ loops-1 ](t1,t2,D5) + dalpha ) ;
      #endif
    }
    #ifdef FIT_D5
    dacorr = fparams[ DATA -> map[i].p[0] ] * asq * ( loop[0](t1,t2,fparams[4]) + dalpha ) / M_PI ;
    dalpha = ( 1 + acr * asq ) * ( loop[ 0 ](t1,t2,fparams[4]) + dalpha ) / M_PI ;
    #else
    dacorr = fparams[ DATA -> map[i].p[0] ] * asq * ( loop[0](t1,t2,D5) + dalpha ) / M_PI ;
    dalpha = ( 1 + acr * asq ) * ( loop[ 0 ](t1,t2,D5) + dalpha ) / M_PI ;
    #endif

    // alpha_s and it's correction terms
    df[ DATA -> map[i].p[0] ][i] = dalpha ;
    df[ DATA -> map[i].p[3] ][i] = dacorr ;
    
    // rotation-preserving corrections
    df[ DATA -> map[i].p[1] ][i] = asq *  ( Qref - DATA -> x[i] ) / ( t1 - t2 ) ;
    df[ DATA -> map[i].p[2] ][i] = asq * asq * ( Qref * Qref - DATA -> x[i] * DATA -> x[i] ) / ( t1 - t2 ) ;

    // derivative of perturbative free parameter d_5
    //df[ DATA -> map[i].p[4] ][i] = pow( a_pi , 5 ) ;
  }
  return ;
}

void
alpha_D0_d2f( double **d2f , const void *data , const double *fparams )
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
alpha_D0_guesses( double *fparams ,
		  const struct data_info Data ,
		  const struct fit_info Fit )
{
  fparams[0] = 0.3 ;
  fparams[1] = 0.17 ;
  fparams[2] = -0.019 ;
  fparams[3] = -3.0 ;
  return ;
}
