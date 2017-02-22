/**
   @file pp_aa_ww.c
   @brief pp_aa correlation functions with the Wall-Walls as support

   MASS = 0
   PL = 1
   AL = 2
   PW = 3
   AW = 4
 */
#include "gens.h"

enum { PWPL , AWAL , PWAL , AWPL , PWPW , AWAW , PWAW , AWPW } ;

double
fpp_aa_ww( const struct x_desc X , const double *fparams , const size_t Npars )
{
  switch( Npars ) {
  case PWPL : // P^W P^L
    return fparams[1] * fparams[3] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case AWAL : // A^W A^L
    return fparams[2] * fparams[4] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case PWAL : // P^W A^L
    return fparams[3] * fparams[2] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case AWPL : // A^L P^L
    return fparams[1] * fparams[4] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case PWPW : // P^W P^W
    return fparams[3] * fparams[3] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case AWAW : // A^W A^W
    return fparams[4] * fparams[4] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case PWAW : // P^W A^W
    return fparams[3] * fparams[4] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case AWPW : // A^W P^W
    return fparams[4] * fparams[3] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  }
  fprintf( stderr , "Should never get here fppaa %zu\n" , Npars ) ;
  exit(-1) ;
  return sqrt(-1) ;
}

void
pp_aa_ww_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fpp_aa_ww( X , fparams , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
pp_aa_ww_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double bck = DATA -> LT[i] - t ;
    const double fwd = exp( -fparams[ 0 ] * t ) ;
    const double bwd = exp( -fparams[ 0 ] * bck ) ;

    // initialise everything to zero
    for( j = 0 ; j < 5 ; j++ ) df[j][i] = 0.0 ;
    
    switch( DATA -> map[i].bnd ) {
    case PWPL :
      // derivative wrt mass
      df[0][i] = -fparams[1] * fparams[3] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[1][i] = fparams[3] * ( fwd + bwd ) ;
      df[3][i] = fparams[1] * ( fwd + bwd ) ;
      break ;
    case AWAL :
      // derivative wrt mass
      df[0][i] = -fparams[2] * fparams[4] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[2][i] = fparams[4] * ( fwd + bwd ) ;
      df[4][i] = fparams[2] * ( fwd + bwd ) ;
      break ;
    case PWAL :
      // derivative wrt mass
      df[0][i] = -fparams[3] * fparams[2] * ( t * fwd - bck * bwd ) ;
      // derivative wrt amplitudes
      df[2][i] = fparams[3] * ( fwd - bwd ) ;
      df[3][i] = fparams[2] * ( fwd - bwd ) ;
      break ;
    case AWPL :
      // derivative wrt mass
      df[0][i] = -fparams[4] * fparams[1] * ( t * fwd - bck * bwd ) ;
      // derivative wrt amplitudes
      df[1][i] = fparams[4] * ( fwd - bwd ) ;
      df[4][i] = fparams[1] * ( fwd - bwd ) ;
      break ;
    case PWPW :
      // derivative wrt mass
      df[0][i] = -fparams[3] * fparams[3] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[3][i] = 2 * fparams[3] * ( fwd + bwd ) ;
      break ;
    case AWAW :
      // derivative wrt mass
      df[0][i] = -fparams[4] * fparams[4] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[4][i] = 2 * fparams[4] * ( fwd + bwd ) ;
      break ;
    case PWAW :
      // derivative wrt mass
      df[0][i] = -fparams[3] * fparams[4] * ( t * fwd - bck * bwd ) ;
      // derivative wrt amplitudes
      df[3][i] = fparams[4] * ( fwd + bwd ) ;
      df[4][i] = fparams[3] * ( fwd + bwd ) ;
      break ;
    case AWPW :
      // derivative wrt mass
      df[0][i] = -fparams[4] * fparams[3] * ( t * fwd - bck * bwd ) ;
      // derivative wrt amplitudes
      df[3][i] = fparams[4] * ( fwd - bwd ) ;
      df[4][i] = fparams[3] * ( fwd - bwd ) ;
      break ;
    }
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
pp_aa_ww_d2f( double **d2f , const void *data , const double *fparams )
{
    const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double bck = DATA -> LT[i] - t ;
    const double fwd = exp( -fparams[ 0 ] * t ) ;
    const double bwd = exp( -fparams[ 0 ] * bck ) ;

    // initialise everything to zero
    for( j = 0 ; j < 25 ; j++ ) d2f[j][i] = 0.0 ;

    // 
    switch( DATA -> map[i].bnd ) {
    case PWPL :
      // second derivative wrt mass
      d2f[0][i]  = fparams[1] * fparams[3] * ( t * t * fwd + bck * bck * bwd ) ;
      // non-zero derivatives wrt mass and then amplitudes
      d2f[5][i]  = d2f[1][i]  = -fparams[3] * ( t * fwd + bck * bwd ) ;
      d2f[15][i] = d2f[3][i]  = -fparams[1] * ( t * fwd + bck * bwd ) ;
      // non-zero derivatives wrt amplitudes
      d2f[16][i] = d2f[8][i] = fwd + bwd ; // fparams1, fparams3
      break ;
    case AWAL :
      // second derivative wrt mass
      d2f[0][i]  = fparams[2] * fparams[4] * ( t * t * fwd + bck * bck * bwd ) ;
      // non-zero derivatives wrt mass and then amplitudes
      d2f[10][i] = d2f[2][i]  = -fparams[4] * ( t * fwd + bck * bwd ) ;
      d2f[20][i] = d2f[4][i]  = -fparams[2] * ( t * fwd + bck * bwd ) ;
      // non-zero derivatives wrt amplitudes
      d2f[22][i] = d2f[14][i] = fwd + bwd ;
      break ;
    case PWAL :
      // second derivative wrt mass
      d2f[0][i]  = fparams[3] * fparams[2] * ( t * t * fwd - bck * bck * bwd ) ;
      // non-zero derivatives wrt mass and then amplitudes
      d2f[15][i] = d2f[3][i]  = -fparams[2] * ( t * fwd - bck * bwd ) ;
      d2f[10][i] = d2f[2][i]  = -fparams[3] * ( t * fwd - bck * bwd ) ;
      // non-zero derivatives wrt amplitudes
      d2f[17][i] = d2f[13][i] = fwd - bwd ;
      break ;
    case AWPL :
      // second derivative wrt mass
      d2f[0][i]  = fparams[4] * fparams[1] * ( t * t * fwd - bck * bck * bwd ) ;
      // non-zero derivatives wrt mass and then amplitudes
      d2f[20][i] = d2f[4][i]  = -fparams[1] * ( t * fwd - bck * bwd ) ;
      d2f[5][i]  = d2f[1][i]  = -fparams[4] * ( t * fwd - bck * bwd ) ;
      // non-zero derivatives wrt amplitudes
      d2f[21][i] = d2f[9][i] = fwd - bwd ;
      break ;
      // WALL-WALL contributions
    case PWPW :
      // second derivative wrt mass
      d2f[0][i]  = fparams[3] * fparams[3] * ( t * t * fwd + bck * bck * bwd ) ;
      // non-zero derivatives wrt mass and then amplitudes
      d2f[15][i]  = d2f[3][i]  = -2 * fparams[3] * ( t * fwd + bck * bwd ) ;
      // non-zero derivatives wrt amplitudes
      d2f[18][i] = 2 * ( fwd + bwd ) ; // fparams3, fparams3
      break ;
    case AWAW :
      // second derivative wrt mass
      d2f[0][i]  = fparams[4] * fparams[4] * ( t * t * fwd + bck * bck * bwd ) ;
      // non-zero derivatives wrt mass and then amplitudes
      d2f[20][i]  = d2f[4][i]  = -2 * fparams[4] * ( t * fwd + bck * bwd ) ;
      // non-zero derivatives wrt amplitudes
      d2f[24][i] = 2 * ( fwd + bwd ) ; // fparams4, fparams4
      break ;
    case PWAW :
    case AWPW :
      // second derivative wrt mass
      d2f[0][i]  = fparams[3] * fparams[4] * ( t * t * fwd - bck * bck * bwd ) ;
      // non-zero derivatives wrt mass and then amplitudes
      d2f[15][i]  = d2f[3][i]  = -2 * fparams[4] * ( t * fwd - bck * bwd ) ;
      d2f[20][i]  = d2f[4][i]  = -2 * fparams[3] * ( t * fwd - bck * bwd ) ;
      // amplitudes
      d2f[23][i] = d2f[19][i]  = 2 * ( fwd - bwd ) ; // fparams4, fparams4
      break ;
    }
  }
  return ;
}

void
pp_aa_ww_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 0.38 ;
  fparams[1] = 20 ;
  fparams[2] = 3 ;
  fparams[3] = 8000 ;
  fparams[4] = 5000 ;
}
