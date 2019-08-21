/**
   @file nrqcd_exp.c
   @brief fit of NRQCD correlators to determine the kinetic mass
 */
#include "gens.h"

#include "Nder.h"

static double psq[ 25 ] = { 0 } ;

//#define SQRT
//#define TAYLOR
#define P4_AMP
//#define P4_EXP
//#define P4_EXP2
//#define P4_EXPAMP2
//#define P4_EXPAMP

void
set_psq_nrqcd2( const double *p2 ,
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
fnrqcd_exp2( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef SQRT
  return fparams[0] * ( 1. + fparams[3] * psq[ Npars ] )
    * exp( -sqrt( fparams[1]*fparams[1] + psq[ Npars ]*fparams[1]/fparams[2] ) * X.X ) ;
#elif (defined TAYLOR)
  return fparams[0] * ( 1. + fparams[3] * psq[ Npars ] )
    * exp( -( fparams[1] + psq[ Npars ]/(2*fparams[2]) ) * X.X  ) ;
#elif (defined P4_EXP)
  return fparams[0] * ( 1. + fparams[3] * psq[ Npars ] )
    * exp( -( fparams[1] + psq[ Npars ]/(2*fparams[2]) + fparams[4]*psq[ Npars ] * psq[ Npars ] ) * X.X ) ;
#elif (defined P4_AMP)
  return fparams[0] * ( 1. + fparams[3] * psq[ Npars ] + fparams[4]*psq[ Npars ] * psq[ Npars ])
    * exp( -( fparams[1] + psq[ Npars ]/(2*fparams[2]) ) * X.X ) ;
#elif (defined P4_EXPAMP)
  return fparams[0] * ( 1. + fparams[3] * psq[ Npars ] + fparams[4]*psq[ Npars ] * psq[ Npars ])
    * exp( -( fparams[1] + psq[ Npars ]/(2*fparams[2]) + fparams[5]*psq[ Npars ] * psq[ Npars ]) * X.X ) ;
#elif (defined P4_EXP2)
  return fparams[0] * ( 1. + fparams[3] * psq[ Npars ] )
    * exp( -( fparams[1] + psq[ Npars ]/(2*fparams[2]) * ( 1 - psq[ Npars ]/8. ) ) * X.X ) ;
#elif (defined P4_EXPAMP2)
  return fparams[0] * ( 1. + fparams[3] * psq[ Npars ] + fparams[4]*psq[Npars]*psq[Npars] )
    * exp( -( fparams[1] + psq[ Npars ]/(2*fparams[2]) * ( 1 - psq[ Npars ]/8. ) ) * X.X ) ;
#endif
}

void
nrqcd_exp2_f( double *f , const void *data , const double *fparams )
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
    f[i] = fnrqcd_exp2( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
nrqcd_exp2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double t = DATA -> x[i] ;
    const size_t bnd = DATA -> map[i].bnd ;
    
    const double M0 = fparams[ DATA -> map[i].p[1] ] ;
    const double MK = fparams[ DATA -> map[i].p[2] ] ;

#ifdef SQRT
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) ;
    const double root = sqrt( M0*M0 + psq[bnd]*M0/MK ) ;
    const double expfac = exp( -root * t ) ;
    df[DATA -> map[i].p[0]][i] = ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) * expfac ;
    df[DATA -> map[i].p[1]][i] = -A * ( 2*M0 + psq[bnd]/MK ) * t * expfac / ( 2*root ) ;
    df[DATA -> map[i].p[2]][i] =  A * ( psq[bnd]*M0/(MK*MK) ) * t * expfac / ( 2*root ) ;
    df[DATA -> map[i].p[3]][i] =  psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
#elif (defined TAYLOR)
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) ;
    const double root = M0 + psq[bnd]/(2*MK) ;
    const double expfac = exp( -root * t ) ;
    df[DATA -> map[i].p[0]][i] = ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) * expfac ;
    df[DATA -> map[i].p[1]][i] = -A * t * expfac ;
    df[DATA -> map[i].p[2]][i] =  A * psq[bnd] * t * expfac / ( 2*MK*MK ) ;
    df[DATA -> map[i].p[3]][i] =  psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
#elif (defined P4_EXP)
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) ;
    const double root = M0 + psq[bnd]/(2*MK) + psq[bnd]*psq[bnd]*fparams[DATA -> map[i].p[4]] ;
    const double expfac = exp( -root * t ) ;
    df[DATA -> map[i].p[0]][i] =  expfac ;
    df[DATA -> map[i].p[1]][i] = -A * t * expfac ;
    df[DATA -> map[i].p[2]][i] =  A * psq[bnd] * t * expfac / ( 2*MK*MK ) ;
    df[DATA -> map[i].p[3]][i] =  psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
    df[DATA -> map[i].p[4]][i] = -A * t * psq[bnd] * psq[bnd] * expfac ;
#elif (defined P4_AMP)
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] + psq[bnd]*psq[bnd]*fparams[ DATA -> map[i].p[4] ]) ;
    const double root = M0 + psq[bnd]/(2*MK) ;
    const double expfac = exp( -root * t ) ;
    df[DATA -> map[i].p[0]][i] = ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] + psq[bnd]*psq[bnd]*fparams[ DATA -> map[i].p[4] ]) * expfac ;
    df[DATA -> map[i].p[1]][i] = -A * t * expfac ;
    df[DATA -> map[i].p[2]][i] =  A * psq[bnd] * t * expfac / ( 2*MK*MK ) ;
    df[DATA -> map[i].p[3]][i] =  psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
    df[DATA -> map[i].p[4]][i] =  psq[bnd] * psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
#elif (defined P4_EXPAMP)
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] + psq[bnd]*psq[bnd]*fparams[ DATA -> map[i].p[4] ]) ;
    const double root = M0 + psq[bnd]/(2*MK) + psq[bnd]*psq[bnd]*fparams[DATA -> map[i].p[5]] ;
    const double expfac = exp( -root * t ) ;
    df[DATA -> map[i].p[0]][i] = expfac ;
    df[DATA -> map[i].p[1]][i] = -A * t * expfac ;
    df[DATA -> map[i].p[2]][i] =  A * psq[bnd] * t * expfac / ( 2*MK*MK ) ;
    df[DATA -> map[i].p[3]][i] =  psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
    df[DATA -> map[i].p[4]][i] =  psq[bnd] * psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
    df[DATA -> map[i].p[5]][i] = -A * t * psq[bnd] * psq[bnd] * expfac ;
#elif (defined P4_EXP2)
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] ) ;
    const double root = M0 + psq[bnd]/(2*MK) * ( 1 - psq[bnd]/8. ) ;
    const double expfac = exp( -root * t ) ;
    df[DATA -> map[i].p[0]][i] =  expfac ;
    df[DATA -> map[i].p[1]][i] = -A * t * expfac ;
    df[DATA -> map[i].p[2]][i] =  A * psq[bnd] * ( 1-psq[bnd]/8. ) * t * expfac / ( 2*MK*MK ) ;
    df[DATA -> map[i].p[3]][i] =  psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
#elif (defined P4_EXPAMP2)
    const double A  = fparams[ DATA -> map[i].p[0] ] * ( 1 + psq[bnd]*fparams[ DATA -> map[i].p[3] ] + psq[bnd]*psq[bnd]*fparams[ DATA -> map[i].p[4] ]) ;
    const double root = M0 + psq[bnd]/(2*MK) * ( 1 - psq[bnd]/8. ) ;
    const double expfac = exp( -root * t ) ;
    df[DATA -> map[i].p[0]][i] =  expfac ;
    df[DATA -> map[i].p[1]][i] = -A * t * expfac ;
    df[DATA -> map[i].p[2]][i] =  A * psq[bnd] * ( 1-psq[bnd]/8. ) * t * expfac / ( 2*MK*MK ) ;
    df[DATA -> map[i].p[3]][i] =  psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
    df[DATA -> map[i].p[4]][i] =  psq[bnd] * psq[bnd] * fparams[ DATA -> map[i].p[0] ] * expfac ;
#endif
  }

  return ;
}

void
nrqcd_exp2_d2f( double **d2f , const void *data , const double *fparams )
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
nrqcd_exp2_guesses( double *fparams ,
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
