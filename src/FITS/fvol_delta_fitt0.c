/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + D*e^{-\sqrt{X}*L }   
 */
#include "gens.h"

#include "fake.h"

#include "Nder.h"
#include <gsl/gsl_sf_bessel.h>

#define COMPUTE_ZOMEGA

// choices for param2
//#define ASQ_TERM
//#define A_TERM
//#define ASQM_TERM

//#define CUTA

//#define CUTL

// GeV values
//static const double mpi = 0.1348 ;
//static const double mk  = 0.4942 ;
//static const double mpi = 0.1349766 ;
//static const double mk  = 0.497614 ;

static const double mpi = 0.1349768 ;
static const double mk = 0.495011 ;


static const double chiral_f = 0.0924 ; // f in GeV
static const double chiral_mu = 0.77 ;
static const double mOct = 0.8043 ;
static const double mDec = 1.1152 ;

// PDG values
static const double mOmega = 1.6695 ; //1.67245 ; //1.6695 ;
static const double mXi = 0.5*(1.31486+1.32171) ; // average over the charged and neutral
static const double mXistar = 0.5*(1.5316+1.5344) ;

// regensberg value in GeV^-2 used here but not so big of a deal
static const double musq = 8*0.5391144626032651*(chiral_mu*chiral_mu) ; 
static const double PhiNucleonCont = 8*0.5391144626032651*(mOct*mOct) ;
static const double PhiDeltaCont = 8*0.5391144626032651*(mDec*mDec) ;
static const double phif = 8*0.5391144626032651*chiral_f*chiral_f ;

#ifdef CUTA
  #define NENSEMBLES (24)
#elif (defined CUTL)
  #define NENSEMBLES (25)
#else
  #define NENSEMBLES (27)
#endif

// lookup table of bubbles and FV terms
struct LUT {
  double b11[ NENSEMBLES+1 ] ;
  double b21[ NENSEMBLES+1 ] ;
  double b22[ NENSEMBLES+1 ] ;
  double IQ0msq[ NENSEMBLES ] ;
  double IQ2[ NENSEMBLES ] ;
  double IMeta[ NENSEMBLES ] ;
  double D1[ NENSEMBLES+1 ] ;
  double D2[ NENSEMBLES+1 ] ;

  // eta
  double IQ0msq_eta[ NENSEMBLES ] ;
  double IQ2_eta[ NENSEMBLES ] ;

  // pion
  double fvpi0[ NENSEMBLES ] ;
  double fvpi2[ NENSEMBLES ] ;
  
  double IMetaSub ;
} ;

static struct LUT lut ;

// need to give this an error
static double Phi3[ NENSEMBLES+1 ] = { };

struct resampled *phi3_res = NULL ;

// GMOR (4*\phi_K-\phi_\pi)/3.
static const double PhiEta[ NENSEMBLES+1 ] = {
#ifndef CUTA
  0.7803, 1.0473, 1.1778,
#endif
#ifndef CUTL
  0.7496,
#endif
  0.7557, 0.9412, 1.1481, 1.1352, 1.2444,
  0.7559, 0.9612, 1.1273, 1.1244, 1.2648, 1.3741,
#ifndef CUTL
  0.7459,
#endif
  0.7410, 0.9647, 1.1335, 1.2916, 1.3833,
  0.7636, 1.0013, 1.2184, 1.3634,
  0.7439, 0.9926,
  -1
  //8*t0*(4*mk*mk-mpi*mpi)/3
};
static const double PhiCascade[ NENSEMBLES+1 ] = {
#ifndef CUTA
  5.954, 6.401, 6.663,
#endif
#ifndef CUTL 
  6.216,
#endif
  5.879, 6.256, 6.722, 6.661, 6.876,
  6.151, 6.426, 6.814, 6.879, 7.002, 7.184,
#ifndef CUTL
  6.488,
#endif
  5.935, 6.521, 6.889, 7.154, 7.302,
  6.297, 6.814, 7.090, 7.322,
  5.988, 6.675,
  -1
  //8*t0*mXi*mXi
} ;
static const double PhiCasStar[ NENSEMBLES+1 ] = {
#ifndef CUTA
  8.620, 9.317, 9.435,
#endif
#ifndef CUTL
  9.257,
#endif
  8.691, 8.943, 9.352, 9.320, 9.524,
  9.021, 9.302, 9.662, 9.662, 9.713, 9.746,
#ifndef CUTL
  9.309,
#endif
  8.887, 9.261, 9.696, 9.944, 10.020,
  9.237, 10.087, 10.006, 10.599, // eww E300 is clearly bad
  9.097, 9.708,
  -1
  //8*t0*mXistar*mXistar
} ;
static const double PhiOmega[ NENSEMBLES+1 ] = {
#ifndef CUTA
  8.620, 9.704, 10.199,
#endif
#ifndef CUTL
  9.257,
#endif
  8.691, 9.476, 10.323, 10.143, 10.778,
  9.021, 9.896, 10.395, 10.411, 11.073, 11.363,
#ifndef CUTL
  9.309,
#endif
  8.887, 9.812, 10.609, 11.290, 11.624,
  9.237, 10.218, 11.134, 11.811,
  9.097, 10.047,
  -1
  //8*t0*mOmega*mOmega
} ;

static const double MKL[ NENSEMBLES+1 ] = {
#ifndef CUTA
  5.09, 5.45, 7.40,
#endif
#ifndef CUTL
  4.33,
#endif
  5.83, 6.13, 6.47, 9.68, 9.88,
  5.13, 5.45, 7.12, 8.56, 11.77, 11.94,
#ifndef CUTL
  4.33,
#endif
  6.44, 6.91, 7.23, 10.00, 15.29,
  5.07, 5.45, 7.66, 11.88,
  5.22, 5.63,
  -1
} ;

static const double MPIL[ NENSEMBLES+1 ] = {
#ifndef CUTA
  5.09, 4.0, 4.10,
#endif
#ifndef CUTL
  4.33,
#endif
  5.83, 4.93, 3.88, 5.83, 4.59,
  5.13, 4.32, 4.42, 5.31, 5.33, 3.80,
#ifndef CUTL
  4.33,
#endif
  6.44, 5.40, 4.43, 4.16, 4.00,
  5.07, 4.17, 4.14, 4.22,
  5.22, 4.21,
  -1
} ;

static const double MetaL[ NENSEMBLES+1 ] = {
#ifndef CUTA
  5.09, 5.86, 8.21,
#endif
#ifndef CUTL
  4.33,
#endif
  5.73, 6.48, 7.13, 10.65, 11.09,
  5.13, 5.77, 9.39, 7.82, 13.23, 13.61,
#ifndef CUTL
  4.33,
#endif
  6.44, 7.35, 7.95, 11.30, 17.50,
  5.07, 5.82, 8.51, 13.50,
  5.22, 6.03,
  -1
} ;

static const double MCasL[ NENSEMBLES+1 ] = {
#ifndef CUTA
  14.052, 14.568, 19.548,
#endif
#ifndef CUTL
  12.463,
#endif
  16.237, 16.698, 17.254, 25.805, 26.078 ,
  14.662, 14.928, 19.316, 23.179, 31.136, 36, // still a guess for D452
#ifndef CUTL
  12.698,
#endif
  18.235, 19.104, 19.603, 26.590, 40.205 ,
  14.568, 15.173, 20.531, 31.277,
  14.803, 15.629,
  -1
} ;

static const double MCsStL[ NENSEMBLES+1 ] = {
#ifndef CUTA
  16.908, 17.4936, 23.262,
#endif
#ifndef CUTL
  15.201,
#endif
  19.741, 19.965, 20.352, 30.523, 30.730,
  17.757, 18.032, 23.916, 28.699, 36.672, 39.0,
#ifndef CUTL
  15.205,
#endif
  22.316, 22.766, 23.256, 31.347, 47.098,
  17.645, 18.461, 24.390, 37.632,
  18.246, 18.848,
  -1
} ;

static const double MOML[ NENSEMBLES+1 ] = {
#ifndef CUTA
  16.908, 17.8536, 24.1856,
#endif
#ifndef CUTL
  15.2007744,
#endif
  19.7414528, 20.5504, 21.3824, 31.8432, 32.5760256,
  17.7569312, 18.5248, 23.744, 28.5408, 39.1547584, 39.51744,
#ifndef CUTL
  15.2047552,
#endif
  22.3161984, 23.4341136, 24.3264, 33.4016, 50.7264,
  17.6448, 18.5808, 25.728, 39.7248,
  18.2464, 19.1744 ,
  -1
} ;

////////////////////////////// FVOL FUNCYIONS AND THINGS /////////////////////
static inline
double mu( const double z ,
	   const double mRL ,
	   const double mQL ,
	   const double mBL )
{
  return sqrt( z*mBL*mBL + (1-z)*mQL*mQL - (1-z)*z*mBL*mBL ) ;
}

static double
integrand( const double z ,
	   const double mRL ,
	   const double mQL ,
	   const double mBL )
{
  const double r2 = 1.4142135623730951 , r3 = 1.7320508075688772 , r5 = 2.23606797749979 ;
  const double r6 = 2.449489742783178 , r8 = 2.8284271247461903 , r10 = 3.1622776601683795 ;
  const double r11 = 3.3166247903554 , r12 = 3.4641016151377544 , r13 = 3.605551275463989 ;
  const double muz = mu( z , mRL , mQL , mBL ) ; 
  return (  +6*gsl_sf_bessel_K0(    muz  ) 
	   +12*gsl_sf_bessel_K0( r2*muz  )/r2
	    +8*gsl_sf_bessel_K0( r3*muz  )/r3
	    +6*gsl_sf_bessel_K0(  2*muz  )/2
	   +24*gsl_sf_bessel_K0( r5*muz  )/r5
	   +24*gsl_sf_bessel_K0( r6*muz  )/r6
	   +12*gsl_sf_bessel_K0( r8*muz  )/r8
	   +30*gsl_sf_bessel_K0(  3*muz  )/3
	   +24*gsl_sf_bessel_K0( r10*muz )/r10
	   +24*gsl_sf_bessel_K0( r11*muz )/r11
	    +8*gsl_sf_bessel_K0( r12*muz )/r12
	   +24*gsl_sf_bessel_K0( r13*muz )/r13
	    ) ;
}

static double
intz_indep( const double mRL ,
	    const double mQL )
{
  const double r2 = 1.4142135623730951 , r3 = 1.7320508075688772 , r5 = 2.23606797749979 ;
  const double r6 = 2.449489742783178 , r8 = 2.8284271247461903 , r10 = 3.1622776601683795 ;
  const double r11 = 3.3166247903554 , r12 = 3.4641016151377544 , r13 = 3.605551275463989 ;
  return 2*(mQL/(mRL*mRL-mQL*mQL))*(   6*gsl_sf_bessel_K1( mQL     )
				     +12*gsl_sf_bessel_K1( r2*mQL  )/r2
				      +8*gsl_sf_bessel_K1( r3*mQL  )/r3
				      +6*gsl_sf_bessel_K1( 2*mQL   )/2
				     +24*gsl_sf_bessel_K1( r5*mQL  )/r5
				     +24*gsl_sf_bessel_K1( r6*mQL  )/r6
				     +12*gsl_sf_bessel_K1( r8*mQL  )/r8
				     +30*gsl_sf_bessel_K1( 3*mQL   )/3
				     +24*gsl_sf_bessel_K1( r10*mQL )/r10
				     +24*gsl_sf_bessel_K1( r11*mQL )/r11
				      +8*gsl_sf_bessel_K1( r12*mQL )/r12
				     +24*gsl_sf_bessel_K1( r13*mQL )/r13 ) ;
}

typedef struct {
  double a , b , eps , whole, fa , fb , fm ;
  int rec ;
  double mRL , mQL , mBL ;
} simp_args ;

static inline double
inner_simp( const double a ,
	    const double b ,
	    const double fa ,
	    const double fb ,
	    const double fm )
{
  return (b-a)*( fa + 4*fm + fb )/6. ; 
}

// adaptive routine adapted fromm wikipedia
// https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method
static double
simpson( const simp_args Simp )
{  
  register const double m  = (Simp.a+Simp.b)/2 ;
  register const double lm = (Simp.a+m)/2 ;
  register const double rm = (m+Simp.b)/2 ;

  const double flm = integrand( lm , Simp.mRL , Simp.mQL , Simp.mBL ) ;
  const double frm = integrand( rm , Simp.mRL , Simp.mQL , Simp.mBL ) ;
  const double left  = inner_simp( Simp.a/2 , Simp.b/2, Simp.fa , Simp.fm , flm ) ;
  const double right = inner_simp( Simp.a/2 , Simp.b/2, Simp.fm , Simp.fb , frm ) ;

  const double delta = left + right - Simp.whole ;
  
  if( Simp.rec < 0 ) {
    printf( "out of recursions %d %e\n" , Simp.rec , delta ) ;
    return sqrt(-1) ;
  }

  if( fabs( delta ) < 15*Simp.eps ) {
    return left + right + delta/15 ;
  }

  const simp_args SimpLeft  = { .a = Simp.a , .b = m , .eps = Simp.eps/2 ,
				.whole = left , .fa = Simp.fa , .fb = Simp.fm ,
				.fm = flm , .rec = Simp.rec-1 ,
				.mRL = Simp.mRL , .mQL = Simp.mQL , .mBL = Simp.mBL } ;
  const simp_args SimpRight = { .a = m , .b = Simp.b , .eps = Simp.eps/2 ,
				.whole = right , .fa = Simp.fm , .fb = Simp.fb ,
				.fm = frm , .rec = Simp.rec-1 ,
				.mRL = Simp.mRL , .mQL = Simp.mQL , .mBL = Simp.mBL } ;
  
  return simpson( SimpLeft ) +	simpson( SimpRight ) ;
}

// simpsons adaptive routine
static double
fv_IQR( const double mRL ,
	const double mQL ,
	const double mBL )
{
  const double a = 0 , b = 1 , h = b-a ;
  const double fa = integrand( a , mRL , mQL , mBL ) ;
  const double fb = integrand( b , mRL , mQL , mBL ) ;
  const double fm = integrand( (a+b)/2 , mRL , mQL , mBL ) ;  
  const double whole = (h/6)*(fa+4*fm+fb) ;
  const simp_args SimpInitial  = { .a = a , .b = 1 , .eps = 1E-12 ,
				   .whole = whole , .fa = fa , .fb = fb ,
				   .fm = fm , .rec = 35 ,
				   .mRL = mRL , .mQL = mQL , .mBL = mBL } ;
  return 2*( simpson( SimpInitial ) - intz_indep( mRL , mQL )) ;
}

void
set_phi3( const size_t sample_idx , const bool is_avg )
{
  for( size_t i = 0 ; i < NENSEMBLES+1 ; i++ ) {
    if( is_avg ) {
      Phi3[i] = phi3_res[i].avg ;
    } else {
      Phi3[i] = phi3_res[i].resampled[ sample_idx ] ;
    }
  } 
}

static inline double
fv_correction_k1( const double mQL )
{
  const double r2 = 1.4142135623730951 , r3 = 1.7320508075688772 , r5 = 2.23606797749979 ;
  const double r6 = 2.449489742783178 , r8 = 2.8284271247461903 , r10 = 3.1622776601683795 ;
  const double r11 = 3.3166247903554 , r12 = 3.4641016151377544 , r13 = 3.605551275463989 ;
  return ( +6*gsl_sf_bessel_K1( mQL     )
	  +12*gsl_sf_bessel_K1( r2*mQL  )/r2
	   +8*gsl_sf_bessel_K1( r3*mQL  )/r3
	   +6*gsl_sf_bessel_K1( 2*mQL   )/2
	  +24*gsl_sf_bessel_K1( r5*mQL  )/r5
	  +24*gsl_sf_bessel_K1( r6*mQL  )/r6
	  +12*gsl_sf_bessel_K1( r8*mQL  )/r8
	  +30*gsl_sf_bessel_K1( 3*mQL   )/3
	  +24*gsl_sf_bessel_K1( r10*mQL )/r10
	  +24*gsl_sf_bessel_K1( r11*mQL )/r11
	   +8*gsl_sf_bessel_K1( r12*mQL )/r12
	  +24*gsl_sf_bessel_K1( r13*mQL )/r13 )/mQL ;
}

static inline double
fv_correction_k2( const double mQL )
{
  const double r2 = 1.4142135623730951 , r3 = 1.7320508075688772 , r5 = 2.23606797749979 ;
  const double r6 = 2.449489742783178 , r8 = 2.8284271247461903 , r10 = 3.1622776601683795 ;
  const double r11 = 3.3166247903554 , r12 = 3.4641016151377544 , r13 = 3.605551275463989 ;
  return (   6*gsl_sf_bessel_Kn( 2 , mQL     )
	   +12*gsl_sf_bessel_Kn( 2 , r2*mQL  )/2
	    +8*gsl_sf_bessel_Kn( 2 , r3*mQL  )/3
	    +6*gsl_sf_bessel_Kn( 2 , 2*mQL   )/4
	   +24*gsl_sf_bessel_Kn( 2 , r5*mQL  )/5
	   +24*gsl_sf_bessel_Kn( 2 , r6*mQL  )/6
	   +12*gsl_sf_bessel_Kn( 2 , r8*mQL  )/8
	   +30*gsl_sf_bessel_Kn( 2 , 3*mQL   )/9
	   +24*gsl_sf_bessel_Kn( 2 , r10*mQL )/10
	   +24*gsl_sf_bessel_Kn( 2 , r11*mQL )/11
	    +8*gsl_sf_bessel_Kn( 2 , r12*mQL )/12
	   +24*gsl_sf_bessel_Kn( 2 , r13*mQL )/13 )/(mQL*mQL) ;
}

static inline double
IQ0( const double phiQ , const double sub , const double MQL )
{
  register double fv = 0.0 ;
  if( MQL != -1 ) {
    fv += 4*phiQ*fv_correction_k1( MQL ) ;
  } 
  return ( phiQ*log(phiQ/musq) - sub + fv )/(16*M_PI*M_PI) ;
}

static inline double
msqIQ0( const double phi , const double sub , const double MQL )
{
  register double fv = 0.0 ;
  if( MQL != -1.0 ) {
    fv += 4*phi*phi*fv_correction_k1( MQL ) ;
  } 
  return (phi*phi*log(phi/musq) - sub + fv )/(16*M_PI*M_PI) ;
}

static inline double
mqIQ2( const double phi , const double sub , const double MQL )
{
  register double fv = 0.0 ;
  if( MQL != -1.0 ) {
    fv += 4*phi*phi*fv_correction_k2( MQL ) ;
  } 
  return ((phi*phi*log(phi/musq) - sub)/4. - fv )/(16*M_PI*M_PI) ;
}

static inline double
PQRsq( const double phiB , const double phiR, const double phiMQ )
{
  return phiB/4. - (phiR + phiMQ)/2. + pow(phiR-phiMQ,2)/(4.*phiB) ;
}

static inline double
ERsq( const double phiB , const double phiR, const double phiMQ )
{
  return phiR + PQRsq( phiB, phiR, phiMQ ) ;
}

static inline double
IQR( const double phiB , const double phiR , const double phiQ , const double PQR , const double fv , const bool is_cascade )
{
  register double gam = 0.0 ;
  // check that B!=R as the log would complain
  if( is_cascade == true ) {
    gam = -((PhiNucleonCont-PhiDeltaCont)/PhiDeltaCont)*log(fabs(PhiNucleonCont-PhiDeltaCont)/PhiNucleonCont ) ;
  }
  register const double i_pqr = sqrt(fabs(PQR)) ;
  register const double mB = sqrt( phiB ) ;
  register double imfac = -2*(i_pqr/mB)*atan2( 2*i_pqr*mB/(phiQ+phiR) , 1 - phiB/(phiQ+phiR ) ) ;
  
  return ( gam
	   + 0.5*( (phiR-phiQ)/phiB - 1. )*log(phiQ/phiR)
	   + imfac
	   + fv
	   ) /(16.*M_PI*M_PI) ;
}

static inline double
delta1( const double M , const double D )
{
  return -( M*(2*M+D)/pow(M+D,2))*log( D*(2*M+D)/(M*M) ) ;
}

static inline double
delta2( const double M , const double D )
{
  return M/(2*M+D) + M*( 2*M*M + 2*D*M + D*D )/( (2*M+D)*pow(M+D,2) )*log( (D*(2*M+D))/(M*M)) ;
}

static inline double
aQR( const double phiOm , const double phiCas , const double phiNuc , const double phiDel , const double phiMk )
{
  const double mOm  = sqrt( phiOm ) ;
  const double mCas = sqrt( phiCas ) ;
  
  const double M = sqrt( phiNuc ) ;
  const double D = sqrt( phiDel ) - M ;
  register const double d1 = delta1( M , D ) ;
  register const double pd1dD = (delta1( M , D+1E-6 )-delta1( M , D-1E-6 ))/(2E-6) ;
  register const double pd1dN = (delta1( M+1E-6 , D )-delta1( M-1E-6 , D ))/(2E-6) ;
  register const double b1 = pow( 2*M+D , 4 )/(16*M*pow(M+D,3)) ;
  return (b1*D*D*( (mOm - M - D)*(d1+D*pd1dD) - (mCas-M)*(D*pd1dD-D*pd1dN+d1*(M+D)/M) )
	  + D*phiMk*b1*delta2( M , D )
	  )/(16*M_PI*M_PI) ;
}

static double
bubble1( const double phiB , const double phiR , const double phiQ , const double MQL , const double MRL , const double fv )
{
  // up to a constant
  const double Er = sqrt( ERsq( phiB , phiR, phiQ ) ) ;
  const double Psq = PQRsq( phiB , phiR , phiQ ) ;

  const double mR = sqrt(phiR) ;
  const double mB = sqrt(phiB) ;
  const double p1 = pow(mR+mB,2)*(phiR-phiB)/( 24.*pow(phiB,1.5) ) ;
  const double IQ1 = IQ0( phiQ , 0 , MQL ) ;
  const double IQ2 = IQ0( phiR , 0 , MRL ) * phiQ/phiR ;
  const double p3 = ( Er+mR)*Psq*IQR( phiB , phiR , phiQ , Psq , fv , true ) ;

  const double aqr = aQR( phiB , phiR , PhiNucleonCont , PhiDeltaCont , phiQ ) ;
  return (p1*(IQ1-IQ2) - p3/3.+ 2*aqr/3.) ; 
}

// needs the cascade* and the eta
static double
bubble2( const double phiB , const double phiR , const double phiQ , const double MQL , const double MRL , const double fv )
{
  const double Er  = sqrt( ERsq( phiB , phiR, phiQ ) ) ;
  const double Psq = PQRsq( phiB , phiR , phiQ ) ;
  const double IQ1 = IQ0( phiQ , 0 , MQL ) ;
  const double IQ2 = IQ0( phiR , 0 , MRL ) * phiQ/phiR ;

  const double mB = sqrt(phiB) ;
  const double mR = sqrt(phiR) ;
  
  const double p1 = ( phiB + phiR )*( phiR - phiB )*(IQ1-IQ2)/(18.*phiB*mR) ;
  const double p2 = ( phiR*phiR + phiB*phiB + 12*phiR*phiB)*( phiR-phiB )*(IQ1-IQ2)/(36.*pow(phiB,1.5)*phiR) ;
  const double p3 = ( (pow(mB+mR,2)/(9.*phiR)) )*(2*Er*(Er-mR)+5*phiR)*Psq*
    IQR( phiB, phiR, phiQ, Psq , fv , false )/(Er+mR) ;

  return -p1 + p2 -p3 ;
}

static struct LUT
setlutcont( const double t0 )
{
  struct LUT lutcont = lut ;
  {
    const double phiOmega   = 8*t0*mOmega*mOmega ;
    const double phiCascade = 8*t0*mXi*mXi ;
    const double phiCasStar = 8*t0*mXistar*mXistar ;

    const double phi3 = 8*t0*mk*mk ;
    const double phi5 = 8*t0*(4*mk*mk-mpi*mpi)/3 ;
    
    const size_t i = NENSEMBLES ;
    const double eps = 1E-6 ;
    const double mOm = sqrt( phiOmega ) ;
    const double pl = pow( mOm + eps , 2 ) , mn = pow( mOm - eps , 2 ) ;
    lutcont.b11[i] = bubble1( phiOmega , phiCascade , phi3 , -1 , -1 , 0 )/phif ;
    lutcont.b21[i] = bubble2( phiOmega , phiCasStar , phi3 , -1 , -1 , 0 )/phif ;
    lutcont.b22[i] = bubble2( phiOmega , phiOmega   , phi5 , -1 , -1 , 0 )/phif ;
    
    lutcont.IMetaSub = IQ0( phi5 , 0.0 , -1 ) ;

    lutcont.D1[ i ] = sqrt( phiOmega )/phif*\
      ( bubble1( pl , phiCascade , phi3 , -1 , -1 , 0.0 )
	- bubble1( mn , phiCascade , phi3 , -1 , -1 , 0.0 ) ) / (2*eps) ;
    lutcont.D2[ i ] = (1/3.)*sqrt( phiOmega )/phif*\
      ( + bubble2( pl , phiCasStar , phi3 , -1 , -1 , 0.0 )
	- bubble2( mn , phiCasStar , phi3 , -1 , -1 , 0.0 )
	+ bubble2( pl , pl , phi3 , -1 , -1 , 0.0 )
	- bubble2( mn , mn , phi3 , -1 , -1 , 0.0 ) ) / (2*eps) ;
  }
  return lutcont ;
}

void
init_phi3( const size_t Nboots )
{  
  // need to give this an error
  static double y[ NENSEMBLES+1 ] = {
#ifndef CUTA
    0.7803, 0.9045, 0.9550,
#endif
#ifndef CUTL
    0.7496,
#endif
    0.7557, 0.8420, 0.9462, 0.9365, 0.9866,
    0.7559, 0.8556, 0.9355, 0.9353, 0.9999, 1.0366,
#ifndef CUTL
    0.7459,
#endif
    0.7410, 0.8539, 0.9382, 1.0126, 1.0556 ,
    0.7636, 0.8797, 0.9860, 1.0558,
    0.7439, 0.8657 ,
    -1
  };
  static double dy[ NENSEMBLES+1 ] = {
#ifndef CUTA
    0.0070, 0.0071, 0.0047,
#endif
#ifndef CUTL
    0.0051,
#endif
    0.0056, 0.0050, 0.0057, 0.0029, 0.0033,
    0.0048, 0.0039, 0.0019, 0.0017, 0.0020, 0.0017,
#ifndef CUTL
    0.0086,
#endif
    0.0053, 0.0024, 0.0029, 0.0019, 0.0012,
    0.0038, 0.0044, 0.0026, 0.0017,
    0.0051, 0.0044 ,
    0.0
  };
  phi3_res = generate_fake_boot( NENSEMBLES+1, Nboots, y , dy ) ;

  set_phi3( 0 , true ) ;
  for( size_t i = 0 ; i < NENSEMBLES+1 ; i++ ) {
    fprintf( stdout , "PHI3 set --> %zu %e %e\n" , i , phi3_res[i].avg , phi3_res[i].err ) ;
  }

  // precompute all bubbles
  fprintf( stdout , "Precomputing bubbles\n" ) ;
  for( size_t i = 0 ; i < NENSEMBLES ; i++ ) {

    const double eps = 1E-6 ;
    const double mOm = sqrt(PhiOmega[ i ]) ;
    const double pl = pow( mOm + eps , 2 ) , mn = pow( mOm - eps , 2 ) ;
    
    const double fv1 = fv_IQR( MCasL[i]  , MKL[i]   , MOML[i] ) ;
    
    lut.b11[i] = bubble1( PhiOmega[ i ] , PhiCascade[ i ] , Phi3[ i ]   , MKL[ i ]   , MCasL[i]  , fv1 )/phif ;


    const double fv2 = fv_IQR( MCsStL[i] , MKL[i]   , MOML[i] ) ;
    lut.b21[i] = bubble2( PhiOmega[ i ] , PhiCasStar[ i ] , Phi3[ i ]   , MKL[ i ]   , MCsStL[i] , fv2 )/phif ;

    const double fv3 = fv_IQR( MOML[i]   , MetaL[i] , MOML[i] ) ; 
    lut.b22[i] = bubble2( PhiOmega[ i ] , PhiOmega[ i ]   , PhiEta[ i ] , MetaL[ i ] , MOML[i]   , fv3 )/phif ;

    lut.IQ0msq[ i ] = msqIQ0( Phi3[i] , 0 , MKL[i] ) ;
    lut.IQ2[ i ]    = mqIQ2( Phi3[i] , 0 , MKL[i] ) ;
    
    // eta lookup table values
    lut.IMeta[i]        = IQ0( PhiEta[i] , 0.0 , MetaL[i] ) ;
    lut.IQ0msq_eta[ i ] = msqIQ0( PhiEta[i] , 0 , MetaL[i] ) ;
    lut.IQ2_eta[ i ]    = mqIQ2( PhiEta[i] , 0 , MetaL[i] ) ;

    lut.fvpi0[ i ] = 4*fv_correction_k1( MPIL[i] )/(16*M_PI*M_PI) ;
    lut.fvpi2[ i ] = 4*fv_correction_k2( MPIL[i] )/(16*M_PI*M_PI) ;
    
    
    lut.D1[ i ] = sqrt(8*0.5391144626032651*mOmega*mOmega)/(phif)*		\
      ( bubble1( pl , PhiCascade[ i ] , Phi3[ i ] , MKL[ i ] , MCasL[ i ] , fv1 )
	- bubble1( mn , PhiCascade[ i ] , Phi3[ i ] , MKL[ i ] , MCasL[ i ] , fv1 ) ) / (2*eps) ;

    lut.D2[ i ] = (1/3.)*sqrt(8*0.5391144626032651*mOmega*mOmega)/(phif)*			\
      ( + bubble2( pl , PhiCasStar[ i ] , Phi3[ i ] , MKL[ i ] , MCsStL[ i ] , fv2 )
	- bubble2( mn , PhiCasStar[ i ] , Phi3[ i ] , MKL[ i ] , MCsStL[ i ] , fv2 )
	+ bubble2( pl , pl , Phi3[ i ] , MetaL[ i ] , MOML[ i ] , fv3 )
	- bubble2( mn , mn , Phi3[ i ] , MetaL[ i ] , MOML[ i ] , fv3 ) ) / (2*eps) ; 
  }
}

#ifdef COMPUTE_ZOMEGA
static inline double
Zomega( const double *fparams , const struct LUT lt , const int Npars )
{  
  return 1 + fparams[6]*fparams[6]*lt.D1[Npars] + fparams[7]*fparams[7]*(lt.D2[Npars]) ;
}
#endif

double
ffvol_delta( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double t0 = fparams[13] ;
  const struct LUT lutcont = setlutcont( t0 ) ;

  const double phi3cont = 8*t0*(mk*mk) ;
  const double phi5cont = 8*t0*(4*mk*mk-mpi*mpi)/3. ; // GMOR physical eta
  const double phi2cont = 8*t0*(mpi*mpi) ;
  
  const double c1 = phi2cont*phi2cont ;
  const double c3 = phi3cont*phi3cont ;
  const double c6 = phi3cont-phi2cont ;
  const double c5 = phi5cont*phi5cont ;
  const double c2 = (c1*log(phi2cont/musq)) ;
  
  const double c4 = (c3*log(phi3cont/musq))/(16*M_PI*M_PI);

  const double c10 = (c5*log(phi5cont/musq))/(16*M_PI*M_PI) ;
  const double c11 = phi3cont+phi2cont/2. ;

  const double IQMPI = msqIQ0( X.X , c2 , -1 ) + X.X*X.X*lut.fvpi0[ Npars ] ;
  const double IQ2MPI = mqIQ2( X.X , c2 , -1 ) - X.X*X.X*lut.fvpi2[ Npars ] ;
  

  #ifdef COMPUTE_ZOMEGA
  const double Zom = Zomega( fparams , lut     , Npars      ) ;
  const double Zph = Zomega( fparams , lutcont , NENSEMBLES ) ;  
  #else
  const double Zom = 1.0 , Zph = 1.0  ;
  #endif

  return fparams[0]*( 1 
		      +fparams[1]*( -(8/3.)*( (Phi3[Npars]-X.X) - c6 )
				    //+(16/9.)*( (Phi3[Npars]-X.X)*lut.IMeta[Npars] - c6*lutcont.IMetaSub )/phif
				   )
		      #ifdef A_TERM
		      +fparams[2]*( fparams[0]*(Phi3[Npars]+X.X/2.))
                      #elif (defined ASQ_TERM)
		      +fparams[2]*fparams[0]*fparams[0]
                      #elif (defined ASQM_TERM)
		      +fparams[2]*fparams[0]*fparams[0]*( -(8/3.)*( (Phi3[Npars]-X.X) - c6 ) )
		      -4*fparams[2]*fparams[0]*fparams[0]/3.*( (Phi3[Npars]+X.X/2.) - c11 ) 
                      #else // MK^4
		      +fparams[2]*( Phi3[Npars]*Phi3[Npars] - c3 )
		      #endif
		      +fparams[3]*( (Phi3[Npars] - X.X)*PhiEta[Npars] - c6*phi5cont ) 
		      -fparams[4]*(lut.IQ0msq[Npars]-c4)/phif
		      -fparams[5]*(lut.IQ2[Npars]-c4/4.)/phif
		      +fparams[6]*fparams[6]*( +lut.b11[Npars]/Zom - lutcont.b11[NENSEMBLES]/Zph )
		      +fparams[7]*fparams[7]*(1/3.)*( +lut.b21[Npars]/Zom - lutcont.b21[NENSEMBLES]/Zph
						      +lut.b22[Npars]/Zom - lutcont.b22[NENSEMBLES]/Zph )
		      -fparams[8]*IQMPI/phif
		      -fparams[9]*IQ2MPI/phif
		      -fparams[10]*(lut.IQ0msq_eta[Npars]-c10)/phif
		      -fparams[11]*(lut.IQ2_eta[Npars]-c10/4.)/phif
		      -4*(fparams[12]+fparams[1]/3.)*( (Phi3[Npars]+X.X/2.) - c11 ) 
		      ) ;
}

void
fvol_delta_f( double *f , const void *data , const double *fparams )
{
  const double mul = (1.67245/mOmega) ;
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> Npars ] ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    f[i] = ffvol_delta( X , p , i ) - DATA -> y[i]*mul ;
  }
  return ;
}

// derivatives
void
fvol_delta_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;

  for( size_t i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    size_t j ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( ffvol_delta ,
					   X ,
					   i,
					   fparams ,
					   j ,
					   DATA -> Npars ) ;
    }
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvol_delta_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvol_delta_guesses( double *fparams ,
		    const struct data_info Data ,
		    const struct fit_info Fit )
{
  return ;
}
