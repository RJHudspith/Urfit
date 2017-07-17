/**
   @file ppaa.c
   @brief simultaneous fit over pp aa ap pa correlation functions

   @warning data must be in the above order

   P^L P^L
   A^L A^L
   P^L A^L
   A^L P^L
 */
#include "gens.h"

#include "exp.h"

//#define INDIVIDUAL

enum { PLPL , ALAL , PLAL , ALPL } ;

double
fppaa( const struct x_desc X , const double *fparams , const size_t Npars )
{
  switch( Npars%4 ) {
  case PLPL : // P^L P^L
    return fparams[1] * fparams[1] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case ALAL : // A^L A^L
    return fparams[2] * fparams[2] * ( exp( -fparams[0] * X.X ) +
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case PLAL : // P^L A^L
    return fparams[1] * fparams[2] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  case ALPL : // A^L P^L
    return fparams[2] * fparams[1] * ( exp( -fparams[0] * X.X ) -
				       exp( -fparams[0] * ( X.LT - X.X ) ) ) ;
  }
  fprintf( stderr , "Should never get here fppaa %zu\n" , Npars ) ;
  exit(-1) ;
  return sqrt(-1) ;
}

void
ppaa_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fppaa( X , fparams , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
ppaa_df( double **df , const void *data , const double *fparams )
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

    switch( DATA -> map[i].bnd ) {
    case PLPL :
      // derivative wrt mass
      df[0][i] = -fparams[1] * fparams[1] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[1][i] = 2 * fparams[1] * ( fwd + bwd ) ;
      break ;
    case ALAL :
      // derivative wrt mass
      df[0][i] = -fparams[2] * fparams[2] * ( t * fwd + bck * bwd ) ;
      // derivative wrt amplitudes
      df[2][i] = 2 * fparams[2] * ( fwd + bwd ) ;
      break ;
    case PLAL :
      // derivative wrt mass
      df[0][i] = -fparams[1] * fparams[2] * ( t * fwd - bck * bwd ) ;
      // derivative wrt amplitudes
      df[1][i] = fparams[2] * ( fwd - bwd ) ;
      df[2][i] = fparams[1] * ( fwd - bwd ) ;
      break ;
    case ALPL :
      // derivative wrt mass
      df[0][i] = -fparams[2] * fparams[1] * ( t * fwd - bck * bwd ) ;
      // derivative wrt amplitudes
      df[1][i] = fparams[2] * ( fwd - bwd ) ;
      df[2][i] = fparams[1] * ( fwd - bwd ) ;
      break ;
    }
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
ppaa_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
ppaa_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 0.25 ;
  fparams[1] = 7000 ;
  fparams[2] = 2000 ;

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
