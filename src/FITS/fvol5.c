/**
   @file fvol4.c
   @brief finite volume fit y = A*( 1 + B*(X-cont) + C*e^{-\sqrt{X}*L } + D*a^2)
 */
#include "gens.h"

#include "Nder.h"
#include <gsl/gsl_sf_bessel.h>

//#define HIGH_LOOPS

//#define ORDERA

static const double mpi = 0.1348 ;
static const double mK  = 0.4942 ;
static const double t0 = 0.538473314101305 ;

static const double phipicont = 8*t0*mpi*mpi ;
static const double phiKcont  = 8*t0*mK*mK ;

static const double chiral_f = 0.0924 ; // f in GeV
static const double chiral_mu = 0.77 ;
static const double musq = 8*t0*(chiral_mu*chiral_mu) ;
static const double phif = 8*t0*chiral_f*chiral_f ;

static const double phiK[ 16 ] = {
    0.7803, 0.9045, 0.9550,
    0.7496, 0.7557, 0.8420, 0.9462,
    0.7559, 0.9355,
    0.7459, 0.7410, 0.8539, 0.9382,1.0125617911199998,
    0.7636,
    phiKcont
} ;

static const double MKL[ 23 ] = {
  // A653, A654 , B650
  5.08632, 5.45088 , 7.4208,
  // U103, H101, H102, H105, X150
  4.3269264, 5.82944, 6.12608, 6.4736, 8.055268, 9.6663216,
  // B450, X451
  5.14016, 7.1212,
  // H200, N202, N293, N200
  4.3312, 6.44352, 6.91296, 7.23408, 10.0032,
  // N300
  5.07312,
  -1, -1, -1, -1, -1, -1
} ;

static const double asq[ 23 ] = {
  0.2532819643738654,0.2532819643738654,0.2532819643738654,
  0.19152593500567394,0.19152593500567394,0.19152593500567394,0.19152593500567394,0.19152593500567394,0.19152593500567394,
  0.14965075255622193,0.14965075255622193,
  0.10603283349102183,0.10603283349102183,0.10603283349102183,0.10603283349102183,0.10603283349102183,
  0.06370463879342395,

  0.2532819643738654, 0.19152593500567394, 0.14965075255622193, 0.10603283349102183, 0.06370463879342395,
  0.0
  
} ;

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
	   +24*gsl_sf_bessel_K1( r13*mQL )/r13
	   )/mQL ;
}

static inline double
msqIQ0( const double phi , const double sub , const double MQL )
{
  register double fv = 0.0 ;
  if( MQL != -1.0 ) {
#ifdef HIGH_LOOPS
    //fv += 4*phi*phi*fv_correction_k1( MQL ) ;
    fv += 4*fv_correction_k1( MQL ) ;
#else
    //fv += 4*phi*phi*exp( -MQL )*sqrt(M_PI/(2*MQL)) ;
    fv += 4*exp( -MQL ) ; ///MQL ;
#endif
  } 
  return fv/(16*M_PI*M_PI) ; //(phi*phi*log(phi/musq) - sub + fv )/(16*M_PI*M_PI) ;
}

double
ffvol5( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double FV = msqIQ0( phiK[Npars] , log( phiKcont / musq ) , MKL[ Npars ] )/phif ;
 
  return fparams[0] * ( 1 + fparams[1]*(X.X-phipicont)
			+fparams[2]*FV
			#ifdef ORDERA
			+fparams[3]*sqrt(asq[Npars])
			#else
			+fparams[3]*asq[Npars]
			#endif
			) ; 
}

void
fvol5_f( double *f , const void *data , const double *fparams )
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

    f[i] = ffvol5( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvol5_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  for( size_t i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    const size_t *mp = DATA->map[i].p ;

    const double FV = msqIQ0( phiK[i] , log( phiKcont / musq ) , MKL[ i ] )/phif ;
  
    df[ mp[0] ][ i ] = ( 1 + fparams[ mp[1] ]*(X.X-phipicont)
			 +fparams[ mp[2] ]*FV
			 #ifdef ORDERA
			 +fparams[ mp[3] ]*sqrt(asq[i])
			 #else
			 +fparams[ mp[3] ]*asq[i]
			 #endif
			 ) ;
    df[ mp[1] ][ i ] = fparams[mp[0]]*(X.X-phipicont) ;
    df[ mp[2] ][ i ] = fparams[mp[0]]*FV ;
    #ifdef ORDERA
    df[ mp[3] ][ i ] = fparams[mp[0]]*sqrt(asq[i]) ; //fparams[mp[4]]*fparams[mp[4]] ; //asq[i] ;
    #else
    df[ mp[3] ][ i ] = fparams[mp[0]]*asq[i] ; //fparams[mp[4]]*fparams[mp[4]] ; //asq[i] ;
    #endif

    df[ mp[4] ][ i ] = 0.0 ; //2*fparams[mp[0]]*fparams[mp[3]]*fparams[mp[4]] ; //asq[i] ;
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvol5_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvol5_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = -120.0 ;
  fparams[1] = 0.1 ;
  fparams[2] = 1 ;
  fparams[3] = 0 ;
  return ;
}
