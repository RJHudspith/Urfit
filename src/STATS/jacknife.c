/**
   @file jacknife.c
   @brief jacknife resampling
 */
#include "gens.h"

#include "stats.h"
#include "summation.h"

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

      // do the jackknife
      for( k = 0 ; k < N ; k++ ) {
	Input -> Data.x[j].resampled[k] = ( sumx - Input -> Data.x[j].resampled[k] ) * NORM ;
	Input -> Data.y[j].resampled[k] = ( sumy - Input -> Data.y[j].resampled[k] ) * NORM ;
      }

      Input -> Data.x[j].restype = Input -> Data.y[j].restype = JackKnife ;     
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
