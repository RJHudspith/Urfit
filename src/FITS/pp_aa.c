/**
   @file pp_aa.c
   @brief simultaneous fit over pp aa ap pa correlation functions

   @warning data must be in the above order

   P^W P^L
   A^W A^L
   P^W A^L
   A^W P^L
 */
#include "gens.h"

#include "exp.h"

//#define INDIVIDUAL
#define WALL

enum { PWPL , AWAL , PWAL , AWPL } ;

double
fpp_aa( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef WALL
  switch( Npars%4 ) {
  case PWPL : // P^W P^L
    return fparams[1] * fparams[3] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case AWAL : // A^W A^L
    return fparams[2] * fparams[4] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case PWAL : // P^W A^L
    return fparams[3] * fparams[2] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case AWPL : // A^W P^L
    return fparams[1] * fparams[4] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  }
#else
  switch( Npars%4 ) {
  case PWPL : // P^W P^L
    return fparams[1] * fparams[1] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case AWAL : // A^W A^L
    return fparams[2] * fparams[2] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case PWAL : // P^W A^L
  case AWPL : // A^W P^L
    return fparams[1] * fparams[2] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  }
#endif
  fprintf( stderr , "Should never get here fppaa %zu\n" , Npars ) ;
  exit(-1) ;
  return sqrt(-1) ;
}

void
pp_aa_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fpp_aa( X , fparams , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
pp_aa_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double t = DATA -> x[i] ;
    const double bck = DATA -> LT[i] - t ;
    const double fwd = exp( -fparams[ 0 ] * t ) ;
    const double bwd = exp( -fparams[ 0 ] * bck ) ;

    // initialise everything to zero
    for( j = 0 ; j < DATA -> Npars ; j++ ) df[j][i] = 0.0 ;

    #ifdef WALL
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
    }
    #else
    switch( DATA -> map[i].bnd ) {
    case PWPL :
      // derivative wrt mass
      df[0][i] = -fparams[1] * fparams[1] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[1][i] = 2*fparams[1] * ( fwd + bwd ) ;
      break ;
    case AWAL :
      // derivative wrt mass
      df[0][i] = -fparams[2] * fparams[2] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[2][i] = 2*fparams[2] * ( fwd + bwd ) ;
      break ;
    case PWAL :
    case AWPL :
      // derivative wrt mass
      df[0][i] = -fparams[1] * fparams[2] * ( t * fwd - bck * bwd ) ;
      // derivative wrt amplitudes
      df[1][i] = fparams[2] * ( fwd - bwd ) ;
      df[2][i] = fparams[1] * ( fwd - bwd ) ;
      break ;
    }
    #endif
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
pp_aa_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
pp_aa_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 0.38 ;
  fparams[1] = 20 ;
  fparams[2] = 5 ;
  fparams[3] = 10000 ;
  fparams[4] = 2000 ;

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
