/**
   @file pmap.c
   @brief fit function parameter map
 */
#include "gens.h"

struct pmap *
parammap( const struct data_info Data ,
	  const struct fit_info Fit )
{
  // index counters
  size_t i , j , k ;

  // allocate the pmap
  struct pmap *map = malloc( Data.Ntot * sizeof( struct pmap ) ) ;
  size_t n = 0 ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    for( j = 0 ; j < Data.Ndata[i] ; j++ ) {
      map[ n ].p = malloc( Fit.Nparam * sizeof( size_t ) ) ;
      n++ ;
    }
  }

  // set the parameters
  register size_t addition = Fit.Nparam ;
  n = 0 ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    for( k = 0 ; k < Fit.Nparam ; k++ ) {
      if( i == 0 ) {
	for( j = 0 ; j < Data.Ndata[i] ; j++ ) {
	  map[ j + n ].p[k] = k ;
	}
      } else {
	if( Fit.Sims[k] == true ) {
	  for( j = 0 ; j < Data.Ndata[i] ; j++ ) {
	    map[ j + n ].p[k] = k ;
	  }
	} else {
	  for( j = 0 ; j < Data.Ndata[i] ; j++ ) {
	    map[ j + n ].p[k] = addition ;
	  }
	  addition++ ;
	}
      }
    }
    n += Data.Ndata[i] ;
  }
  
  return map ;
}

// free the map structure
void
free_pmap( struct pmap *map ,
	   const size_t n )
{
  size_t i ;
  for( i = 0 ; i < n ; i++ ) {
    free( map[ i ].p ) ;
  }
  free( map ) ;
  return ;
}
