/**
   @file histogram.c
   @brief creates a histogram of the data
 */
#include "gens.h"
#include "gen_ders.h"

#include "resampled_ops.h"
#include "rng.h"
#include "stats.h"

// quicksort comparator
static int
cmp( const void *a ,
     const void *b )
{
  return *(double*)a < *(double*)b ? -1 : \
    ( *(double*)a > *(double*)b ? 1 : 0 ) ;
}

// creates a histogram
struct histogram*
histogram( const double *data ,
	   const size_t Ndata ,
	   const size_t Nbins )
{
  struct histogram *histogram = malloc( Nbins * sizeof( struct histogram ) ) ;

  double *sdata = malloc( Ndata * sizeof( double ) ) ;
  memcpy( sdata , data , Ndata * sizeof( double ) ) ;
  
  // sort the data
  qsort( sdata , Ndata , sizeof( double ) , cmp ) ;

  size_t i , j = 0 ;
  const double binwidth = ( sdata[Ndata-1] - sdata[0] ) / (double)Nbins ;

  // add data in each bin -> binmax is the upper limit of the bin
  double binmax = sdata[0] + binwidth ;
    
  for( i = 0 ; i < Nbins ; i++ ) {
    // set the width
    histogram[i].width = binwidth ;
    histogram[i].binmid = ( i + 0.5 ) * binwidth + sdata[0] ;
    
    histogram[i].count = 0 ;
    while( sdata[j] < binmax && j < Ndata ) {
      j++ ;
      histogram[i].count++ ;
    }

    #ifdef VERBOSE
    printf( "[HIST] %zu %f %zu \n" , i ,
	    histogram[i].binmid , histogram[i].count ) ;
    #endif
    
    binmax += binwidth ;
  }

  free( sdata ) ;

  return histogram ;
}

double
hist_min( size_t *min ,
	  const struct histogram *hist ,
	  const size_t Nbin )
{
  // compute the general finite-difference derivatives
  double y[ Nbin ] , x[ Nbin ] ;
  size_t k ;
  
  const size_t order = 12 > Nbin ? Nbin - 3 : 12 ;
  
  for( k = 0 ; k < Nbin ; k++ ) {
    y[k] = hist[k].count ;
    x[k] = hist[k].binmid ;
  }
      
  double **g = get_ders( y , x , Nbin , order ) ;

#ifdef VERBOSE
  for( k = 0 ; k < Nbin ; k++ ) {
    printf( "Ders %zu :: %f %f %f %f \n" , k , hist[k].binmid , g[k][0] , g[k][1] , g[k][2] ) ;
  }
#endif

  // find the maxima
  size_t max_idx[ Nbin ] , idx = 0 ;
  for( k = 0 ; k < Nbin ; k++ ) {
    max_idx[k] = 0 ;
  }
  
  for( k = 1 ; k < Nbin-1 ; k++ ) {
    // check for a change in derivative around k
    if( ( g[k][1] > 0 && g[k+1][1] < 0 ) ) {
      size_t kmax = k ;
      // g[k+1] is our possible max
      if( g[k][0] < g[k+1][0] ) {
	kmax = k+1 ;
      }
      // make sure that our max has negative curvature
      if( g[kmax][2] < 0 ) {
	max_idx[idx] = kmax ;
	#ifdef VERBOSE
	printf( "Max :: %zu \n" , kmax ) ;
	#endif
	idx++ ;
      }
      // 
    }
  }

  #ifdef VERBOSE
  printf( "idx :: %zu \n" , idx ) ;
  #endif

  double xmin = UNINIT_FLAG ;

  // if there is a single minimum how do we find the middle, we don't
  // we just return one of the edges
  if( idx < 2 ) {
    if( max_idx[0] > Nbin/2 ) {
      *min = 0 ;
      return hist[0].binmid ;
    } else {
      *min = Nbin-1 ;
      return hist[Nbin-1].binmid ;
    }
  }

  // find a minimum point that is bounded by the two maxima
  size_t max , minbest = 123456789 ;
  for( max = 0 ; max < idx-1 ; max++ ) {
    for( k = max_idx[max] ; k < max_idx[max+1] ; k++ ) {
      if( g[k][0] < minbest ) {
	minbest = g[k][0] ;
	*min = k ;
	#ifdef VERBOSE
	printf( "[MIN] bounded :: minbest :: %zu | min :: %zu \n" ,
		minbest , *min ) ;
	#endif
      }
    }
  }

  if( min == 0 ) {
    goto memfree ;
  }

  xmin = ( hist[*min].binmid - g[*min][1] * ( hist[*min].binmid - hist[*min+1].binmid ) /
	   ( g[*min][1] - g[*min+1][1] ) ) ;

 memfree :
  
  for( k = 0 ; k < Nbin ; k++ ) {
    free( g[k] ) ;
  }
  free( g ) ;
      
  return xmin ;
}
