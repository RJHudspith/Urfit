/**
   @file nrqcd_exp.c
   @brief fit of NRQCD correlators to determine the kinetic mass
 */
#include "gens.h"
#include "Nder.h"

static double psq[ 25 ] = { 0 } ;

void
set_psq( const double *p2 ,
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
fnrqcd_exp( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[0] * exp( -( psq[ Npars ]/(2*fparams[1]) ) * X.X ) ;
}

void
nrqcd_exp_f( double *f , const void *data , const double *fparams )
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
    f[i] = fnrqcd_exp( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
nrqcd_exp_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double t = DATA -> x[i] ;
    const size_t bnd = DATA -> map[i].bnd ;
    const size_t p0 = DATA -> map[i].p[0] ;
    const size_t p1 = DATA -> map[i].p[1] ;
    const double expfac = exp( -( psq[bnd]/(2*fparams[p1]) ) * t ) ; 
    df[p0][i] = expfac ;
    df[p1][i] = fparams[p0] *
      ( psq[bnd] / ( 2 * fparams[p1] * fparams[p1] ) * t ) * expfac ;
  }
  return ;
}

void
nrqcd_exp_d2f( double **d2f , const void *data , const double *fparams )
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
nrqcd_exp_guesses( double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit )
{
  size_t i ;
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    fparams[i] = 1. ;
  }
  fparams[1] = 5. ;
}
