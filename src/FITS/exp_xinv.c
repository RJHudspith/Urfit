/**
   exponential Model f = ( A * exp(-lambda * i) - y_i )
 */
#include "gens.h"

#include "fit_chooser.h"
#include "GLS.h"
#include "pade_laplace.h"

double
fexp_xinv( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[0] * ( 1 + fparams[2]/X.X )* exp( -fparams[1] * X.X ) ;
}

void
exp_xinv_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ 3 ] ;
    p[ 0 ] = fparams[ DATA -> map[ i ].p[ 0 ] ] ;
    p[ 1 ] = fparams[ DATA -> map[ i ].p[ 1 ] ] ;
    p[ 2 ] = fparams[ DATA -> map[ i ].p[ 2 ] ] ;
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fexp_xinv( X , p , DATA -> N * 2 ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
exp_xinv_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {    
    const double t = DATA -> x[i] ;
    const double e = exp( -fparams[ DATA-> map[i].p[1] ] * t ) ;

    df[ DATA-> map[i].p[0] ][i] =
      (1+ fparams[ DATA-> map[i].p[2] ]/t)*e ;
    df[ DATA-> map[i].p[1] ][i] = \
      -t * fparams[ DATA-> map[i].p[0] ] *\
      ( 1 + fparams[ DATA-> map[i].p[2] ]/t ) * e ;
    df[ DATA-> map[i].p[2] ][i] = fparams[ DATA-> map[i].p[0] ]*e/t ;
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
exp_xinv_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// guesses using the data, we don't need to be too accurate
void
exp_xinv_guesses( double *fparams ,
	     const struct data_info Data ,
	     const struct fit_info Fit )
{
  // usual counters, p0 and p1 count how many times the fit parameter
  // is used which we average over
  double f[ 2 * Fit.N ] ;
  size_t i , j , shift = 0 , p[ 2*Fit.N ] ;
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    fparams[i] = 0.0 ;
  }
  for( i = 0 ; i < 2*Fit.N ; i++ ) {
    p[i] = 0 ; f[i] = 0.0 ;
  }
  
  // loop each individual data set and fit using least squares to
  // log(y),x
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    size_t N = 0 ;
    for( j = shift ; j < shift + Data.Ndata[i] ; j++ ) {
      if( Data.x[j].avg > Data.LT[j]/2 ) continue ;
      N++ ;
    }

    // set the data
    double x[ N ] , y[ N ] ;
    N = 0 ;
    for( j = shift ; j < shift + Data.Ndata[i] ; j++ ) {
      //if( Data.x[j].avg > Data.LT[j]/2 ) continue ;
      x[ N ] = Data.x[j].avg ;
      y[ N ] = Data.y[j].avg ;
      N++ ;
    }

    #ifdef VERBOSE
    for( j = 0 ; j < N ; j++ ) {
      printf( "%f %f \n" , x[j] , y[j] ) ;
    }
    #endif

    // pade laplace is bad here because we are already in the plateau region!!!!
    if( pade_laplace( f , x , y , N , Fit.N , 3 ) == FAILURE ) {

      fprintf( stderr , "Jamie : do something about this!\n" ) ;
      break ;
    } else {

      for( j = 0 ; j < 2*Fit.N ; j+=2 ) {
	fparams[ Fit.map[shift].p[ j + 0 ] ] +=  f[ j + 0 ] ;
	fparams[ Fit.map[shift].p[ j + 1 ] ] += -f[ j + 1 ] ;
	if( Fit.map[shift].p[ j + 0 ] == j + 0 ) p[j + 0]++ ;
	if( Fit.map[shift].p[ j + 1 ] == j + 1 ) p[j + 1]++ ;
      }

    }

    shift += Data.Ndata[i] ;
  }

  // normalise the guesses if we have summed multiple ones
  for( j = 0 ; j < 2*Fit.N ; j++ ) {
    if( p[j] > 1 ) {
      fparams[ j ] /= p[j] ;
    }
  }

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
