/**
   @file fsol.c
   @brief fit of finite t correlators


   shared m_0 and A

   different cutoff values and c^2

   Data order is 00 stw1,stw2,nstw1,nstw2
   
 */
#include "gens.h"

#include "Nder.h"

static double psq[ 25 ] = { 0 } ;

// put a P4 term in the exponent
//#define EXP_P4

void
set_psq_sol2( const double *p2 ,
	      const size_t NP2 )
{
  size_t i ;
  for( i = 0 ; i < NP2 ; i++ ) {
    psq[i] = p2[i] ;
    printf( "P2 set %zu %f \n" , i , psq[i] ) ; 
  }
  return ;
}

double
fsol2( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef EXP_P4
  if( Npars > 2 ) {
    const double fwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*(fparams[5] + psq[Npars]*fparams[7])) * X.X ) ;
    const double bwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*(fparams[5] + psq[Npars]*fparams[7])) *( X.LT - X.X ) ) ;
    return fparams[0] * ( 1. + ( fparams[6] ) * psq[ Npars ] )*( fwd + bwd ) ;
  } else {
    const double fwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*(fparams[2] + psq[Npars]*fparams[4])) * X.X ) ;
    const double bwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*(fparams[2] + psq[Npars]*fparams[4])) *( X.LT - X.X ) ) ;
    return fparams[0] * ( 1. + ( fparams[3] ) * psq[ Npars ] )* ( fwd + bwd ) ;
  }
#else
  if( Npars > 2 ) {
    const double fwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[5] ) * X.X ) ;
    const double bwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[5] ) *( X.LT - X.X ) ) ;
    return fparams[0] * ( 1. + ( fparams[6] + fparams[7] * psq[ Npars] ) * psq[ Npars ] )*( fwd + bwd ) ;
  } else {
    const double fwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[2] ) * X.X ) ;
    const double bwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[2] ) *( X.LT - X.X ) ) ;
    return fparams[0] * ( 1. + ( fparams[3] + fparams[4] * psq[ Npars] ) * psq[ Npars ] )* ( fwd + bwd ) ;
  }
#endif
}

void
sol2_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fsol2( X , fparams , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
sol2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double t = DATA -> x[i] ;
    const size_t bnd = DATA -> map[i].bnd ;

    const double M0 = fparams[ 1 ] ;
    
    for( int j = 0 ; j < 8 ; j++ ) {
      df[j][i] = 0.0 ;
    }

#ifdef EXP_P4
    if( bnd > 2 ) {
      const double MK = fparams[ 5 ] ;
      const double root = sqrt( M0*M0 + psq[bnd]*(MK + psq[bnd]*fparams[7])) ;
      const double fwd = exp( -root * t ) ;
      const double bwd = +exp( -root * ( DATA -> LT[i] - t )) ;
      const double A  = fparams[0] * ( 1 + psq[bnd]*( fparams[6] ) ) ;
      df[0][i] = ( 1 + psq[bnd]*(fparams[6])) * ( fwd + bwd ) ;
      df[1][i] = -A * M0 * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / root ;
      df[5][i] = -A * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
      df[6][i] = psq[bnd] * fparams[0] * (fwd+bwd) ;
      df[7][i] = -A * psq[bnd] * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
    } else {
      const double MK = fparams[ 2 ] ;
      const double root = sqrt( M0*M0 + psq[bnd]*(MK + psq[bnd]*fparams[4])) ;
      const double fwd = exp( -root * t ) ;
      const double bwd = +exp( -root * ( DATA -> LT[i] - t )) ;
      const double A  = fparams[0] * ( 1 + psq[bnd]*( fparams[3] ) ) ;
      df[0][i] = ( 1 + psq[bnd]*(fparams[3]) ) * ( fwd + bwd ) ;
      df[1][i] = -A * M0 * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / root ;
      df[2][i] = -A * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
      df[3][i] = psq[bnd] * fparams[0] * (fwd+bwd) ;
      df[4][i] = -A * psq[bnd] * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
    }
#else
    if( bnd > 2 ) {
      const double MK = fparams[ 5 ] ;
      const double root = sqrt( M0*M0 + psq[bnd]*MK ) ;
      const double fwd = exp( -root * t ) ;
      const double bwd = +exp( -root * ( DATA -> LT[i] - t )) ;
      const double A  = fparams[0] * ( 1 + psq[bnd]*( fparams[6] + psq[bnd]*fparams[7] ) ) ;
      df[0][i] = ( 1 + psq[bnd]*(fparams[6] + psq[bnd]*fparams[7]) ) * ( fwd + bwd ) ;
      df[1][i] = -A * M0 * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / root ;
      df[5][i] = -A * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
      df[6][i] = psq[bnd] * fparams[0] * (fwd+bwd) ;
      df[7][i] = psq[bnd] * psq[bnd] * fparams[0] * (fwd+bwd) ;
    } else {
      const double MK = fparams[ 2 ] ;
      const double root = sqrt( M0*M0 + psq[bnd]*MK ) ;
      const double fwd = exp( -root * t ) ;
      const double bwd = +exp( -root * ( DATA -> LT[i] - t )) ;
      const double A  = fparams[0] * ( 1 + psq[bnd]*( fparams[3] + psq[bnd]*fparams[4] ) ) ;
      df[0][i] = ( 1 + psq[bnd]*(fparams[3] + psq[bnd]*fparams[4]) ) * ( fwd + bwd ) ;
      df[1][i] = -A * M0 * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / root ;
      df[2][i] = -A * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
      df[3][i] = psq[bnd] * fparams[0] * (fwd+bwd) ;
      df[4][i] = psq[bnd] * psq[bnd] * fparams[0] * (fwd+bwd) ;
    }
#endif
  }

  return ;
}

void
sol2_d2f( double **d2f , const void *data , const double *fparams )
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
sol2_guesses( double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit )
{
  size_t i ;
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    fparams[i] = 1E5 ;
  }
  fparams[1] = 0.5 ;
  fparams[2] = 9. ;
}
