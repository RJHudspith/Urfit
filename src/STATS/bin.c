/**
   @file bin.c
   @brief bin our raw data
 */
#include "gens.h"
#include "stats.h"

// do some binning
static int
bin_sample( struct resampled *Data ,
	    const size_t Binning )
{
  const size_t NSAMPLES_NEW = Data -> NSAMPLES / Binning ;

  fprintf( stdout , "NSAMPLES_NEW %zu\n" , NSAMPLES_NEW ) ;
  
  double *bin = malloc( NSAMPLES_NEW * sizeof( double ) ) , ave = 0.0 ;
  size_t idx , data_idx = 0 , k ;
  for( k = 0 ; k < NSAMPLES_NEW ; k++ ) {
    // reset local binning index
    idx = 0 ;
    ave = 0.0 ;
    // loop binned values
    while( idx < Binning ) {
      ave += Data -> resampled[ data_idx ] ;
      data_idx++ ;
      idx++ ;
    }
    bin[k] = ave / (double)Binning ;
  }

  // reallocate and copy over
  Data -> resampled =							\
    realloc( Data -> resampled , NSAMPLES_NEW * sizeof( double ) ) ;
  Data -> NSAMPLES = NSAMPLES_NEW ;

  for( k = 0 ; k < NSAMPLES_NEW ; k++ ) {
    Data -> resampled[k] = bin[k] ;
  }

  free( bin ) ;
  return SUCCESS ;
}

// perform binning on the data
int
bin_data( struct input_params *Input )
{
  size_t i , j = 0 , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    if( Input -> Traj[i].Bin == 1 ) {
      continue ;
    }
    //
    if( Input -> Data.x[j].NSAMPLES < Input -> Traj[i].Bin ||
	Input -> Data.y[j].NSAMPLES < Input -> Traj[i].Bin ) {
      fprintf( stderr , "[STATS] Binning factor greater than NSAMPLES %zu > %zu\n" , Input -> Traj[i].Bin , Input -> Data.y[j].NSAMPLES ) ;
      return FAILURE ;
    }
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      if( Input -> Data.x[j].restype != Raw ||
	  Input -> Data.y[j].restype != Raw ) {
	fprintf( stderr , "[STATS] I won't bin non-raw data\n" ) ;
	return FAILURE ;
      }
    }
    j += shift ;
  }

  // do the binning
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    if( Input -> Traj[i].Bin < 2 ) {
      shift += Input -> Data.Ndata[i] ;
      continue ;
    }
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      // bin both x and y individually
      bin_sample( &Input -> Data.x[j] , Input -> Traj[i].Bin ) ;
      bin_sample( &Input -> Data.y[j] , Input -> Traj[i].Bin ) ;

      compute_err( &Input -> Data.x[j] ) ;
      compute_err( &Input -> Data.y[j] ) ;

      //#ifdef VERBOSE
      printf( "BINNING :: %f %f %f %f \n" ,
	      Input -> Data.x[j].avg , Input -> Data.x[j].err ,
	      Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
      //#endif
    }
    shift = j ;
  }

  return SUCCESS ;
}
