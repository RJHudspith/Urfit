/**
   @file fsol.c
   @brief fit of charmionium correlators to determine the kinetic mass
 */
#include "gens.h"

#include "Nder.h"

static double psq[ 25 ] = { 0 } ;

//#define SINH
//#define EXP
#define P4

void
set_psq_sol( const double *p2 ,
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
fsol( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double fwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[2] ) * X.X ) ;
#ifdef EXP
  const double bwd = 0.0 ;
#elif (defined SINH)
  const double bwd = -exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[2] ) *( X.LT - X.X ) ) ;
#else
  const double bwd = exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[2] ) *( X.LT - X.X ) ) ;
#endif
  
#ifdef P4
  return fparams[0] * ( 1. + ( fparams[3] + fparams[4] * psq[ Npars] ) * psq[ Npars ] )
    * ( fwd + bwd ) ;
#else
  return fparams[0] * ( 1. + ( fparams[3] ) * psq[ Npars ] )
    * ( fwd + bwd ) ;
#endif
}

void
sol_f( double *f , const void *data , const double *fparams )
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
    f[i] = fsol( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
sol_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double t = DATA -> x[i] ;
    const size_t bnd = DATA -> map[i].bnd ;
    
    const double M0 = fparams[ DATA -> map[i].p[1] ] ;
    const double MK = fparams[ DATA -> map[i].p[2] ] ;

    const double root = sqrt( M0*M0 + psq[bnd]*MK ) ;
    const double fwd = exp( -root * t ) ;
#ifdef EXP
    const double bwd = 0.0 ;
#elif (defined SINH)
    const double bwd = -exp( -root * ( DATA -> LT[i] - t )) ;
#else // cosh
    const double bwd = +exp( -root * ( DATA -> LT[i] - t )) ;
#endif
    
#ifdef P4
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*( fparams[ DATA -> map[i].p[3] ] + psq[bnd]*fparams[ DATA -> map[i].p[4] ] ) ) ;
    df[DATA -> map[i].p[0]][i] = ( 1 + psq[bnd]*(fparams[ DATA -> map[i].p[3] ] + psq[bnd]*fparams[ DATA -> map[i].p[4] ]) ) * ( fwd + bwd ) ;
    df[DATA -> map[i].p[1]][i] = -A * M0 * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / root ;
    df[DATA -> map[i].p[2]][i] = -A * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
    df[DATA -> map[i].p[3]][i] = psq[bnd] * fparams[ DATA -> map[i].p[0] ] * (fwd+bwd) ;
    df[DATA -> map[i].p[4]][i] = psq[bnd] * psq[bnd] * fparams[ DATA -> map[i].p[0] ] * (fwd+bwd) ;
#else
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) ;
    df[DATA -> map[i].p[0]][i] = ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) * ( fwd + bwd ) ;
    df[DATA -> map[i].p[1]][i] = -A * M0 * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / root ;
    df[DATA -> map[i].p[2]][i] = -A * psq[bnd] * ( t * fwd + ( DATA -> LT[i] - t ) * bwd ) / (2*root) ;
    df[DATA -> map[i].p[3]][i] = psq[bnd] * fparams[ DATA -> map[i].p[0] ] * (fwd+bwd) ;
#endif
  }

  return ;
}

void
sol_d2f( double **d2f , const void *data , const double *fparams )
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
sol_guesses( double *fparams ,
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
