/**
   @file jacknife.c
   @brief jacknife resampling

   Is actually a Bias-corrected double super jacknife or something
 */
#include "gens.h"

#include "raw.h"
#include "summation.h"

// do we want to do a second-level bias correction? I think we do
// it seems like it does little to the final result when I tested
// against a correlator
// Taken from "Double jackknife bias-corrected estimators" by berg - 1991
#define DJBCORRECT

// the jackknife error definition incorporates the code for the average and variance
void
jackknife_error( struct resampled *replicas )
{
  // average and variance come from the usual defs
  double ave , err ;
  average( &ave , &err , replicas -> resampled , replicas -> NSAMPLES ) ;

  err = sqrt( err * ( 1.0 - 1.0 / (double)replicas -> NSAMPLES ) ) ;
  replicas -> avg = ave ;
  replicas -> err_hi = ave + err ;
  replicas -> err_lo = ave - err ;
  replicas -> err    = err ;
  return ;
}

// perform a bootstrap resampling on the data
void
jackknife_full( struct input_params *Input )
{
  size_t i , j = 0 , k , shift = 0 ;

  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      const size_t Ntot = Input -> Data.x[j].NSAMPLES <= Input -> Data.Nboots ?
	Input -> Data.Nboots : Input -> Data.x[j].NSAMPLES ;
      const size_t Nspill = Input -> Data.x[j].NSAMPLES ;
      
      const size_t N = Input -> Data.x[j].NSAMPLES ;
      const double NORM = 1.0 / ( N - 1.0 ) ;

      double *x = malloc( N * sizeof( double ) ) ;
      double *y = malloc( N * sizeof( double ) ) ;
      
      for( k = 0 ; k < N ; k++ ) {
	x[k] = Input -> Data.x[j].resampled[k] ;
	y[k] = Input -> Data.y[j].resampled[k] ;
      }

      const double sumx = kahan_summation( x , N ) ;
      const double sumy = kahan_summation( y , N ) ;

      // do the jackknife, perhaps with a bias correction step
      double Bb[ N ] , Bc[ N ] ;
      for( k = 0 ; k < N ; k++ ) {
	Bb[k] = ( sumx - Input -> Data.x[j].resampled[k] ) * NORM ;
	Bc[k] = ( sumy - Input -> Data.y[j].resampled[k] ) * NORM ;
        #ifdef DJBCORRECT
	size_t l ;
	register double sumb = 0.0 , sumc = 0.0 ;
	for( l = 0 ; l < N-1 ; l++ ) {
	  size_t ll = l ;
	  if( l >= k ) ll = l + 1 ;
	  sumb -= ( sumx -
		    Input -> Data.x[j].resampled[k] -
		    Input -> Data.x[j].resampled[ll] )  ;
	  sumc -= ( sumy -
		    Input -> Data.y[j].resampled[k] -
		    Input -> Data.y[j].resampled[ll] ) ;
	}
        Bb[k] = Bb[k]*(N-1.0) + sumb * NORM ;
        Bc[k] = Bc[k]*(N-1.0) + sumc * NORM ;
	#endif
      }
      
      for( k = 0 ; k < N ; k++ ) {
	Input -> Data.x[j].resampled[k] = Bb[k] ;
	Input -> Data.y[j].resampled[k] = Bc[k] ;
      }

      Input -> Data.x[j].restype = Input -> Data.y[j].restype = JackKnife ;     
      jackknife_error( &(Input -> Data.x[j]) ) ;
      jackknife_error( &(Input -> Data.y[j]) ) ;

      #ifdef VERBOSE
      printf( "JACKNIFE %f %f || %f %f \n" ,
	      Input -> Data.x[j].avg , Input -> Data.x[j].err ,
	      Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
      #endif

      // reallocate to a bigger value and put the average at the end
      Input -> Data.x[j].NSAMPLES = Input -> Data.y[j].NSAMPLES = Ntot ;
      Input -> Data.x[j].resampled =		\
	realloc( Input -> Data.x[j].resampled , Ntot * sizeof( double ) ) ;
      Input -> Data.y[j].resampled = \
	realloc( Input -> Data.y[j].resampled , Ntot * sizeof( double ) ) ;
      for( k = N ; k < Ntot ; k++ ) {
	Input -> Data.x[j].resampled[k] = Input -> Data.x[j].avg ;
	Input -> Data.y[j].resampled[k] = Input -> Data.y[j].avg ;
      }
      jackknife_error( &(Input -> Data.x[j]) ) ;
      jackknife_error( &(Input -> Data.y[j]) ) ;
      
      #ifdef VERBOSE
      printf( "JACKNIFE %f %f || %f %f \n" ,
	      Input -> Data.x[j].avg , Input -> Data.x[j].err ,
	      Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
      #endif

      free( x ) ; free( y ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  return ;
}
