/**
   @file cruel_runnings.c
   @brief running code for alpha_s using embedded cash-carp

   This version is thread safe
 */
#include "gens.h"

#include "resampled_ops.h"
#include "stats.h"

// defined in the on shell scheme ...
#define MC (1.6) //1.29 // charm mass  OS Scheme 
#define MB (4.7) //4.19 // bottom mass OS Scheme
#define MT (172.9) // Top mass
#define MZ (91.1876) // z-mass
#define HALLEY_TOL (1.E-14) // tolerance for Halley's method

// LUT for beta function
const static double beta[ 7 ][ 5 ] =
  {
    { 8.753521870054244e-01, 6.459225457199034e-01, 7.198643271535230e-01, 1.172686527092228e+00 , 1.7141380966000839 } , // nf = 0
    { 8.223005393081260e-01, 5.657099420030526e-01, 5.819927502680805e-01, 9.104347415045818e-01 , 1.17538974969005 } , // nf = 1
    { 7.692488916108274e-01, 4.854973382862019e-01, 4.501870001712895e-01, 6.810330551205890e-01 , 0.7442815692062783 } , // nf = 2
    { 7.161972439135290e-01, 4.052847345693511e-01, 3.244470768631501e-01, 4.848422163521834e-01 , 0.416174720228734 } , // nf = 3
    { 6.631455962162307e-01, 3.250721308525004e-01, 2.047729803436621e-01, 3.222229736112991e-01 , 0.18628925534643193 } , // nf = 4
    { 6.100939485189321e-01, 2.448595271356496e-01, 9.116471061282574e-02, 1.935360753098698e-01 , 0.0497041146574356 } , // nf = 5
    { 5.570423008216338e-01, 1.646469234187989e-01, -1.637773232935912e-02, 9.914226985982916e-02 , 0.001357125768857487 }  // nf = 6
  } ;

// embedded cash-karp factors
const static double b21 = 0.2 ;
const static double b31 = 0.075 ;
const static double b32 = 0.225 ;
const static double b41 = 0.3 ;
const static double b42 = -0.9 ;
const static double b43 = 1.2 ;
const static double b51 = -0.20370370370370369 ;
const static double b52 = 2.5 ;
const static double b53 = -2.5925925925925926 ;
const static double b54 = 1.2962962962962963 ;
const static double b61 = 0.029495804398148147 ;
const static double b62 = 0.341796875 ;
const static double b63 = 0.041594328703703706 ;
const static double b64 = 0.40034541377314814 ;
const static double b65 = 0.061767578125 ;
const static double c1 = 0.097883597883597878 ;
const static double c3 = 0.40257648953301128 ;
const static double c4 = 0.21043771043771045 ;
const static double c6 = 0.28910220214568039 ;
const static double dc1 = -0.0042937748015873106 ;
const static double dc3 = 0.018668586093857853 ;
const static double dc4 = -0.034155026830808066 ;
const static double dc5 = -0.019321986607142856 ;
const static double dc6 = 0.039102202145680387 ;

// OK cool, write the adaptive bit
const static double ADAPTIVE_EPS = 1E-12 ;
// Standard shrink and factor from NRC
const static double ADAPTIVE_SHRINK = -0.25 ;
// Standard growth and factor from NRC
const static double ADAPTIVE_GROWTH = -0.20 ;
// define adaptive safe
const static double ADAPTIVE_SAFE = 0.9 ;

// adaptive error conserving
#define ADAPTIVE_ERRCON pow( 5./ADAPTIVE_SAFE , 1./ADAPTIVE_GROWTH )

// fmin and fmax functions for the adaptive
static double adaptfmax( const double a , const double b ) { return ( b < a ? a : b ) ; }
static double adaptfmin( const double a , const double b ) { return ( a < b ? a : b ) ; }

// all hail glorious beta function!
static inline double
beta_function( const double alpha , const size_t nf , const size_t loops )
{
  // try with gotos
  if( loops > 1 ) goto loop2 ;
  return -alpha * alpha * beta[nf][0] ;
 loop2 :
  if( loops > 2 ) goto loop3 ;
  return -alpha * alpha * ( beta[nf][0] + alpha * beta[nf][1] ) ;
 loop3 : 
  if( loops > 3 ) goto loop4 ;
  return -alpha * alpha * ( beta[nf][0] + alpha * ( beta[nf][1] + alpha * beta[nf][2] ) ) ;
 loop4 :
  if( loops > 4 ) goto loop5 ;
  return -alpha * alpha * ( beta[nf][0] + alpha * ( beta[nf][1] + alpha * ( beta[nf][2] + alpha * beta[nf][3] ) ) ) ;
 loop5 :
  return -alpha * alpha * ( beta[nf][0] + alpha * ( beta[nf][1] + alpha * ( beta[nf][2] + alpha * ( beta[nf][3] + alpha * beta[nf][4] ) ) ) ) ;
}

// run with the RK4
static double
RK4_step( const double mu , double *alpha , 
	  const size_t nf , const size_t loops ,
	  const double dlq , const double exp_step )
{
  const double F1 = beta_function( *alpha , nf , loops ) ;
  const double F2 = beta_function( *alpha + 0.5 * dlq * F1 , nf , loops ) ;
  const double F3 = beta_function( *alpha + 0.5 * dlq * F2 , nf , loops ) ;
  const double F4 = beta_function( *alpha + dlq * F3 , nf , loops ) ; 
  *alpha += dlq * ( F1 + 2.0 * F2 + 2.0 * F3 + F4 ) * 0.16666666666666666 ;
  return mu * exp_step ;
}

// this is the cash-karp embedded adaptive runge-kutta step
static void
RK4_adapt( const double mu , double *alpha , 
	   const size_t nf , const size_t loops ,
	   const double dlq , double *err )
{
  const double K1 = beta_function( *alpha , nf , loops ) ;
  const double K2 = beta_function( *alpha + dlq*b21*K1 , nf , loops ) ;
  const double K3 = beta_function( *alpha + dlq*(b31*K1 + b32*K2) , nf , loops ) ;
  const double K4 = beta_function( *alpha + dlq*(b41*K1 + b42*K2 + b43*K3) , nf , loops ) ;
  const double K5 = beta_function( *alpha + dlq*(b51*K1 + b52*K2 + b53*K3 + b54*K4) , nf , loops ) ;
  const double K6 = beta_function( *alpha + dlq*(b61*K1 + b62*K2 + b63*K3 + b64*K4 + b65*K5) , nf , loops ) ;
  // compute the new alpha
  *alpha += dlq*( c1*K1 + c3*K3 + c4*K4 + c6*K6 ) ; 
  *err = dlq*( dc1*K1 + dc3*K3 + dc4*K4 + dc5*K5 + dc6*K6 ) ;
  return ;
}

// driver for the adaptive algorithm
static double
adaptive( const double mu , double *alpha ,
	  const size_t nf , const size_t loops , 
	  double *dlq )
{
  double errmax = 10.0 , a1 ;
  while( errmax > 1.0 ) {
    // initialise our attempts
    a1 = *alpha ;
    RK4_adapt( mu , &a1 , nf , loops , *dlq , &errmax ) ;

    errmax = fabs(errmax) / ADAPTIVE_EPS ;
    if( errmax < 1.0 ) { break ; } 

    // set up a tolerance s.t del_temp doesn't go too crazy, only really used when starting guess is bad
    register const double del_temp =  ADAPTIVE_SAFE * (*dlq) * pow( errmax , ADAPTIVE_SHRINK ) ; 
    register const double tol = 0.1 * (*dlq) ; 
    *dlq = ( 0.0 < del_temp ? adaptfmax( del_temp , tol ) : adaptfmin( del_temp , tol ) ) ;
  }
  // and set the variables to their new values
  *alpha = a1 ;
  const double step = mu * exp( 0.5 * (*dlq) ) ;
  *dlq = errmax > ADAPTIVE_ERRCON ? ADAPTIVE_SAFE * (*dlq) * pow( errmax , ADAPTIVE_GROWTH ) : ADAPTIVE_SAFE * 5.0 * (*dlq) ;
  return step ;
}

// I invert this numerically opposed to using the inverse series so I have numerical forwards
//  -> backwards symmetry
static double
match_up_OS( const double Alpha_Ms ,
	     const double Mh ,
	     const size_t nf ,
	     const size_t loops )
{
  if( loops <= 2 )
    return Alpha_Ms ; // nice result, threshold effects only kick in at two loops
                      // but we always match the running to one lower order in the matching ...
  double h2_term = 0. , h3_term = 0. , h4_term = 0. ;
  // switch for the scheme we are matching to
  if( loops > 2 ) { h3_term = -0.02955201190 ; }
  if( loops > 3 ) { h4_term = -0.1717036285 + ( nf - 1 ) * 0.008465086429 ; }

  //set up guess and new_guess
  double guess = Alpha_Ms , correction = 1.0 ;
  
  // don't need to do so many as halley's method converges like x^3
  size_t counter = 0 ;
  while( fabs( correction ) > HALLEY_TOL ) {
    //function and its first and second derivatives
    const double f_x = guess * ( 1. + guess * ( h2_term + guess * ( h3_term + guess * ( h4_term ) ) ) ) - Alpha_Ms ; 
    const double f_prime = 1. + guess * ( 2. * h2_term + guess * ( 3. * h3_term + guess * ( 4. * h4_term ) ) ) ; 
    const double f_pprime = 2. * h2_term + guess * ( 6. * h3_term + guess * ( 12. * h4_term ) ) ;
    //halley's method	  
    correction = -2. * f_x * f_prime / (2. * f_prime * f_prime - f_x * f_pprime) ;

    //printf( "CORRECTION %d %e \n" , counter , correction ) ;

    guess += correction ;

    // make sure it doesn't continue forever
    if( counter > 50 ) {
      printf( "Match up OS not converging ... Leaving \n" ) ;
      return -1 ;
    }
    counter ++ ;
  }
  return guess ;
} 

// matching down from nf to nf-1 flavours at scale Mh
static inline double
match_down_OS( const double Alpha_Ms ,
	       const double Mh ,
	       const size_t nf ,
	       const size_t loops )
{
  if( loops <= 2 )
    return Alpha_Ms ; // nice result, threshold effects only kick in at two loops
                      // but we always match the running to one lower order in the matching ...
  double h3_term , h4_term ;
  if( loops > 2 ) goto loop3 ;
  return Alpha_Ms ;
 loop3 :
  h3_term = -0.02955201190 ;
  if( loops > 3 ) goto loop4 ;
  return Alpha_Ms * ( 1.0 + Alpha_Ms * Alpha_Ms * h3_term ) ; 
 loop4 : 
  h4_term = -0.1717036285 + ( nf - 1 ) * 0.008465086429 ;
  return Alpha_Ms * ( 1.0 + Alpha_Ms * Alpha_Ms * ( h3_term + Alpha_Ms * h4_term ) ) ; 
}

// run from mu to mu-prime
double
RUN( double mu ,
     const double alpha_mu ,
     const double muprime , 
     const size_t nf ,
     const size_t loops )
{
  // these get initialised if we need them, because the integration is constant
  // in the step size "mu" we need only compute the exponentials here.
  // However, for the adaptive this is not the case.
  double alpha_mup = alpha_mu ;
  int flag = 0 ;

  // perfectly sensible starting point
  double dlq = 0.01 ;
  while( mu < muprime ) {
    mu = adaptive( mu , &alpha_mup , nf , loops , &dlq ) ;
    flag = 1 ;
  }
  // run backwards ...
  double mdlq = -0.01 ;
  while( mu > muprime && flag == 0 ) {
    mu = adaptive( mu , &alpha_mup , nf , loops , &mdlq ) ;
  }

  // final RK4 step to go to exactly the scale we wish to be
  register const double murat = muprime/mu ;
  const double diff = 2.0 * log( murat ) ; 
  mu = RK4_step( mu , &alpha_mup , nf , loops , diff , murat ) ;

  return alpha_mup ;
}

// running from nf=3 to MZ matching at the charm
// I want to automate this for many nfs
double 
run_nf3_2MZ( const double alpha ,
	     const double mu ,
	     const size_t loops )
{
  double alpha_mup = RUN( mu , alpha , MC , 3 , loops ) ;
  alpha_mup = match_up_OS( alpha_mup , MC , 4 , loops ) ;
  alpha_mup = RUN( MC , alpha_mup , MB , 4 , loops ) ;
  alpha_mup = match_up_OS( alpha_mup , MC , 5 , loops ) ;
  return RUN( MB , alpha_mup , MZ , 5 , loops ) ;
}

// running from nf=5 at MZ matching at the charm and back to some value in Nf=3
// I want to automate this for many nfs
double 
run_MZ_2nf3( const double alpha ,
	     const double mu ,
	     const size_t loops )
{
  double alpha_mup = RUN( MZ , alpha , MB , 5 , loops ) ;
  alpha_mup = match_down_OS( alpha_mup , MB , 5 , loops ) ;
  alpha_mup = RUN( MB , alpha_mup , MC , 4 , loops ) ;
  alpha_mup = match_down_OS( alpha_mup , MC , 4 , loops ) ;
  return RUN( MC , alpha_mup , mu , 3 , loops ) ;
}

// run an nf=3 distribution to MZ
struct resampled
run_distribution_nf3_2MZ( struct resampled alpha ,
			  const double mu ,
			  const size_t loops )
{
  struct resampled alpha_amz = init_dist( NULL , alpha.NSAMPLES ,
					  alpha.restype ) ;
  
  size_t i ;
  for( i = 0 ; i < alpha.NSAMPLES ; i++ ) {
    alpha_amz.resampled[i] = run_nf3_2MZ( alpha.resampled[i] , mu , loops ) ;
  }
  alpha_amz.avg = run_nf3_2MZ( alpha.avg , mu , loops ) ;
  compute_err( &alpha_amz ) ;
  
  return alpha_amz ;
}

// test the running of the coupling from running down to nf=3 @ 2 GeV and back
void
test_running( void )
{
  const double amz = 0.1184 ;

  double alpha = run_MZ_2nf3( amz , 2.0 , 4 ) ;

  printf( "%f -> %f -> %f \n" , amz , alpha , run_nf3_2MZ( alpha , 2.0 , 4 ) ) ;
  
  return ;
}
