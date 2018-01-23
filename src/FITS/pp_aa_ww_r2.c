/**
   @file pp_aa_ww.c
   @brief pp_aa correlation functions with the Wall-Walls and arbitrary r^2 as support

   expects 4 WL first in the order (PP,AA,AP,PA)

   then 4 WW second in the order (PP,AA,AP,PA)

   and finally an arbitrary multiple of 4 of r^2 extended ops
 */
#include "gens.h"

double
fpp_aa_ww_r2( const struct x_desc X , const double *fparams , const size_t Npars )
{
  switch( Npars%4 ) {
  case 0 :
    return fparams[ 1 + Npars/2 ] * fparams[3] * ( exp( -fparams[0] * X.X ) +
						     exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case 1 :
    return fparams[ 2 + Npars/2 ] * fparams[4] * ( exp( -fparams[0] * X.X ) +
						     exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case 2 :
    return fparams[ 1 + Npars/2 ] * fparams[3] * ( exp( -fparams[0] * X.X ) -
						     exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case 3 :
    return fparams[ Npars/2 ] * fparams[4] * ( exp( -fparams[0] * X.X ) -
						     exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  }
  fprintf( stderr , "Should never get here fpp_aa_ww_r2 %zu\n" , Npars ) ;
  exit(-1) ;
  return sqrt(-1) ;
}

void
pp_aa_ww_r2_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fpp_aa_ww_r2( X , fparams , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
pp_aa_ww_r2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double bck = DATA -> LT[i] - t ;
    const double fwd = exp( -fparams[ 0 ] * t ) ;
    const double bwd = exp( -fparams[ 0 ] * bck ) ;
    const size_t Npars = DATA -> map[i].bnd ;

    // initialise everything to zero
    for( j = 0 ; j < 5 ; j++ ) df[j][i] = 0.0 ;

    // WW is special
    if( ( Npars > 3 ) && ( Npars < 8 ) ) {
      switch( Npars ) {
      case 4 :
	// derivative wrt mass
	df[0][i] = -fparams[3] * fparams[3] * ( t * fwd + bck * bwd ) ;
	// derivative wrt amplitudes
	df[3][i] = 2 * fparams[3] * ( fwd + bwd ) ;
	break ;
      case 5 :
	// derivative wrt mass
	df[0][i] = -fparams[4] * fparams[4] * ( t * fwd + bck * bwd ) ;
	// derivative wrt amplitudes
	df[4][i] = 2 * fparams[4] * ( fwd + bwd ) ;
	break ;
      case 6 :
	// derivative wrt mass
	df[0][i] = -fparams[4] * fparams[3] * ( t * fwd - bck * bwd ) ;
	// derivative wrt amplitudes
	df[3][i] = fparams[4] * ( fwd + bwd ) ;
	df[4][i] = fparams[3] * ( fwd + bwd ) ;
	break ;
      case 7 :
	// derivative wrt mass
	df[0][i] = -fparams[3] * fparams[4] * ( t * fwd - bck * bwd ) ;
	// derivative wrt amplitudes
	df[3][i] = fparams[4] * ( fwd - bwd ) ;
	df[4][i] = fparams[3] * ( fwd - bwd ) ;
	break ;
      }
    } else {
      switch( Npars%4 ) {
      case 0 :
	df[0][i] = -fparams[1 + Npars/2] * fparams[3] * ( t * fwd + bck * bwd ) ;
	// derivative wrt amplitudes
	df[1 + Npars/2][i] = fparams[3] * ( fwd + bwd ) ;
	df[3][i] = fparams[1 + Npars/2] * ( fwd + bwd ) ;
	break ;
      case 1 :
	df[0][i] = -fparams[2 + Npars/2] * fparams[4] * ( t * fwd + bck * bwd ) ;
	// derivative wrt amplitudes
	df[2 + Npars/2][i] = fparams[4] * ( fwd + bwd ) ;
	df[4][i] = fparams[2 + Npars/2] * ( fwd + bwd ) ;
	break ;
      case 2 :
	df[0][i] = -fparams[1 + Npars/2] * fparams[3] * ( t * fwd - bck * bwd ) ;
	// derivative wrt amplitudes
	df[1 + Npars/2][i] = fparams[3] * ( fwd - bwd ) ;
	df[3][i] = fparams[1 + Npars/2] * ( fwd - bwd ) ;
	break ;
      case 3 :
	df[0][i] = -fparams[Npars/2] * fparams[4] * ( t * fwd - bck * bwd ) ;
	// derivative wrt amplitudes
	df[Npars/2][i] = fparams[4] * ( fwd - bwd ) ;
	df[4][i] = fparams[Npars/2] * ( fwd - bwd ) ;
	break ;
      }
    }
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
pp_aa_ww_r2_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
pp_aa_ww_r2_guesses( double *fparams ,
		     const struct data_info Data ,
		     const struct fit_info Fit )
{
  fparams[0] = 1.0 ;

  size_t i ;
  // tell us about the guesses always use the prior as a guess
  printf( "\n" ) ;
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    if( Fit.Prior[i].Initialised == true ) {
      fparams[i] = Fit.Prior[i].Val ;
    } 
    printf( "[GUESS] Fit param guess %zu -> %f \n" , i , fparams[i] ) ; 
  }
  printf( "\n" ) ;

  return ;
}
