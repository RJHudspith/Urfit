/**
   @file sort.c
   @brief sort our distribution data by x.avg
 */
#include "gens.h"

#include "resampled_ops.h"

//#define VERBOSE

// median of three values
static double
median_of_three( const double A1 ,
		 const double A2 ,
		 const double A3 )
{
  if( A1 < A2 ) {
    if( A2 < A3 ) {
      return A2 ;
      // A3 < A2
    } else {
      if( A1 < A3 ) {
	return A3 ;
      } else {
	return A1 ;
      }
    }
  } else { // A1 >= A2
    if( A1 < A3 ) {
      return A1 ;
    } else { // A1 is biggest
      if( A2 < A3 ) {
	return A3 ;
      } else {
	return A2 ;
      }
    }
  }
}

// swap distributions
static void
swap( struct input_params *Input ,
      struct resampled tempx ,
      struct resampled tempy ,
      const size_t idx1 ,
      const size_t idx2 )
{
  equate( &tempx , Input -> Data.x[idx1] ) ;
  equate( &Input -> Data.x[idx1] , Input -> Data.x[idx2] ) ;
  equate( &Input -> Data.x[idx2] , tempx ) ;

  equate( &tempy , Input -> Data.y[idx1] ) ;
  equate( &Input -> Data.y[idx1] , Input -> Data.y[idx2] ) ;
  equate( &Input -> Data.y[idx2] , tempy ) ;
  
  return ;
}

// perform a quicksort between specific ranges lo and hi
// using 3-way partitioning
static void
quicksort( struct input_params *Input ,
	   struct resampled tempx ,
	   struct resampled tempy ,
	   const size_t lo ,
	   const size_t hi )
{
  if( hi <= lo ) return ;
  
  register size_t i = lo , lt = lo , gt = hi ;
  
  const double v = median_of_three( Input -> Data.x[lo].avg ,
				    Input -> Data.x[(lo+hi)/2].avg ,
				    Input -> Data.x[hi].avg ) ;
  
  while( i <= gt ) {
    if( Input -> Data.x[i].avg < v ) {
      swap( Input , tempx , tempy , lt++ , i++ ) ;
    } else if( Input -> Data.x[i].avg > v ) {
      swap( Input , tempx , tempy , i , gt-- ) ;
    } else {
      i++ ;
    }
  }

  if( (lt) < 1 ) lt = gt ;

  quicksort( Input , tempx , tempy , lo , lt-1 ) ;
  quicksort( Input , tempx , tempy , gt+1 , hi ) ;
  
  return ;
}

// bog-standard quicksort
int
quick_sort_data( struct input_params *Input )
{
  size_t i , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    struct resampled tempx = init_dist( &Input -> Data.x[ shift ] ,
					Input -> Data.x[ shift ].NSAMPLES ,
					Input -> Data.x[ shift ].restype ) ;
    struct resampled tempy = init_dist( &Input -> Data.y[ shift ] ,
					Input -> Data.y[ shift ].NSAMPLES ,
					Input -> Data.y[ shift ].restype ) ;

    quicksort( Input , tempx , tempy ,
	       shift , shift+Input->Data.Ndata[i]-1 ) ;
    shift += Input->Data.Ndata[i] ;
  }

  #ifdef VERBOSE
  // test
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t j ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      fprintf( stdout , "TEST -> %f %f %f %f \n" ,
	       Input -> Data.x[j].avg ,
	       Input -> Data.x[j].err_lo ,
	       Input -> Data.x[j].err_hi ,
	       Input -> Data.y[j].avg ) ;
    }
    shift = j ;
  }
  #endif
  
  return SUCCESS ;
}

// slow insertion sort
int
insertion_sort_data( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
    
  // sort in ascending order the data by x-avg using an insertion sort
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    struct resampled tempx = init_dist( &Input -> Data.x[ shift ] ,
					Input -> Data.x[ shift ].NSAMPLES ,
					Input -> Data.x[ shift ].restype ) ;
    struct resampled tempy = init_dist( &Input -> Data.y[ shift ] ,
					Input -> Data.y[ shift ].NSAMPLES ,
					Input -> Data.y[ shift ].restype ) ;
    
    // insertion sort
    for( j = shift + 1 ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      double x = Input -> Data.x[ j ].avg ;
      equate( &tempx , Input -> Data.x[j] ) ;
      equate( &tempy , Input -> Data.y[j] ) ;
      int hole = (int)j - 1 ;
      while( hole >= (int)shift && Input -> Data.x[hole].avg > x ) {	
	// copy data
	equate( &Input -> Data.x[hole+1] , Input -> Data.x[hole] ) ;
	equate( &Input -> Data.y[hole+1] , Input -> Data.y[hole] ) ;
	hole-- ;
      }
      equate( &Input -> Data.x[hole+1] , tempx ) ;
      equate( &Input -> Data.y[hole+1] , tempy ) ;
      #ifdef VERBOSE
      fprintf( stdout , "Sorting :: %zu %d %zu \n" , j , hole , Input -> Data.Ndata[i] ) ;
      #endif
    }
    shift = j ;    
    free( tempx.resampled ) ;
    free( tempy.resampled ) ;
  }

#ifdef VERBOSE
  // test
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      fprintf( stdout , "TEST -> %f %f %f %f \n" ,
	       Input -> Data.x[j].avg ,
	       Input -> Data.x[j].err_lo ,
	       Input -> Data.x[j].err_hi ,
	       Input -> Data.y[j].avg ) ;
    }
    shift = j ;
  }
#endif
  
  return SUCCESS ;
}
