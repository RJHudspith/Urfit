/**
   @file ppaa.c
   @brief simultaneous fit over pp aa ap pa correlation functions

   @warning data must be in the above order

   P^L P^L
   A^L A^L
   P^L A^L
 */
#include "gens.h"

#include "Nder.h"

//#define PPAP_ONLY

#ifdef PPAP_ONLY
enum { PLPL , PLAL } ;
#else
enum { PLPL , ALAL , PLAL  } ;
#endif

double
fppaa( const struct x_desc X , const double *fparams , const size_t Npars )
{
  double sum = 0 ;
  size_t i ;
#ifdef PPAP_ONLY
  switch( Npars%2 ) {
  case PLPL : // P^L P^L
    for( i = 0 ; i < 3 * X.N ; i+=3 ) {
      sum += fparams[i+1] * fparams[i+1] * ( exp( -fparams[i] * X.X ) +
					     exp( -fparams[i] * ( X.LT - X.X ) ) ) ;
    }
    return sum ;
  case PLAL : // P^L A^L
    for( i = 0 ; i < 3 * X.N ; i+=3 ) {
      sum += fparams[i+1] * fparams[i+2] * ( exp( -fparams[i] * X.X ) -
					     exp( -fparams[i] * ( X.LT - X.X ) ) ) ;
    }
    return sum ;
  }
#else
  switch( Npars%3 ) {
  case PLPL : // P^L P^L
    for( i = 0 ; i < 3 * X.N ; i+=3 ) {
      sum += fparams[i+1] * fparams[i+1] * ( exp( -fparams[i] * X.X ) +
					     exp( -fparams[i] * ( X.LT - X.X ) ) ) ;
    }
    return sum ;
  case ALAL : // A^L A^L
    for( i = 0 ; i < 3 * X.N ; i+=3 ) {
      sum += fparams[i+2] * fparams[i+2] * ( exp( -fparams[i] * X.X ) +
					     exp( -fparams[i] * ( X.LT - X.X ) ) ) ;
    }
    return sum ;
  case PLAL : // P^L A^L
    for( i = 0 ; i < 3 * X.N ; i+=3 ) {
      sum += fparams[i+1] * fparams[i+2] * ( exp( -fparams[i] * X.X ) -
					     exp( -fparams[i] * ( X.LT - X.X ) ) ) ;
    }
    return sum ;
  }
#endif
  fprintf( stderr , "Should never get here fppaa %zu\n" , Npars ) ;
  exit(-1) ;
  return sqrt(-1) ;
}

void
ppaa_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> N * 3 ] ;
    for( j = 0 ; j < DATA -> N * 3 ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
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
 
    // initialise everything to zero
    for( j = 0 ; j < DATA -> Npars ; j++ ) df[j][i] = 0.0 ;

    #ifdef PPAP_ONLY
    for( j = 0 ; j < 3*DATA -> N ; j+=3 ) {
      const double t = DATA -> x[i] ;
      const double bck = DATA -> LT[i] - t ;
      const double fwd = exp( -fparams[ j ] * t ) ;
      const double bwd = exp( -fparams[ j ] * bck ) ;
      
      switch( DATA -> map[i].bnd ) {
      case PLPL :
	// derivative wrt mass
	df[j+0][i] = -fparams[j+1] * fparams[j+1] * ( t * fwd + bck * bwd ) ;
	// derivative wrt amplitudes
	df[j+1][i] = 2 * fparams[j+1] * ( fwd + bwd ) ;
	break ;
      case PLAL :
	// derivative wrt mass
	df[j+0][i] = -fparams[j+1] * fparams[j+2] * ( t * fwd - bck * bwd ) ;
	// derivative wrt amplitudes
	df[j+1][i] = fparams[j+2] * ( fwd - bwd ) ;
	df[j+2][i] = fparams[j+1] * ( fwd - bwd ) ;
	break ;
      }
    }
    #else
    for( j = 0 ; j < 3*DATA -> N ; j+=3 ) {
      const double t = DATA -> x[i] ;
      const double bck = DATA -> LT[i] - t ;
      const double fwd = exp( -fparams[ j ] * t ) ;
      const double bwd = exp( -fparams[ j ] * bck ) ;
      
      switch( DATA -> map[i].bnd ) {
      case PLPL :
	// derivative wrt mass
	df[j+0][i] = -fparams[j+1] * fparams[j+1] * ( t * fwd + bck * bwd ) ;
	// derivative wrt amplitudes
	df[j+1][i] = 2 * fparams[j+1] * ( fwd + bwd ) ;
	break ;
      case ALAL :
	// derivative wrt mass
	df[j+0][i] = -fparams[j+2] * fparams[j+2] * ( t * fwd + bck * bwd ) ;
	// derivative wrt amplitudes
	df[j+2][i] = 2 * fparams[j+2] * ( fwd + bwd ) ;
	break ;
      case PLAL :
	// derivative wrt mass
	df[j+0][i] = -fparams[j+1] * fparams[j+2] * ( t * fwd - bck * bwd ) ;
	// derivative wrt amplitudes
	df[j+1][i] = fparams[j+2] * ( fwd - bwd ) ;
	df[j+2][i] = fparams[j+1] * ( fwd - bwd ) ;
	break ;
      }
    }
    #endif
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
  fparams[1] = 36 ;
  fparams[2] = 4 ;

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
