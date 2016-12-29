/**
   @file pmap.c
   @brief fit function parameter map
 */
#include "gens.h"

// set the parameter map
struct pmap *
parammap( const size_t Nparams ,
	  const size_t Nsims ,
	  const size_t *Ndata ,
	  const bool *sim_params )
{
  // compute the full map size
  size_t i , j , k , n = 0 ;
  for( i = 0 ; i < Nsims ; i++ ) {
    n += Ndata[i] ;
  }

  // allocate the pmap
  struct pmap *map = malloc( n * sizeof( struct pmap ) ) ;
  n = 0 ;
  for( i = 0 ; i < Nsims ; i++ ) {
    for( j = 0 ; j < Ndata[i] ; j++ ) {
      map[ n ].p = malloc( Nparams * sizeof( double ) ) ;
      n++ ;
    }
  }

  // set the parameters
  register size_t addition = Nparams ;
  n = 0 ;
  for( i = 0 ; i < Nsims ; i++ ) {
    for( k = 0 ; k < Nparams ; k++ ) {
      //
      if( i == 0 ) {
	for( j = 0 ; j < Ndata[i] ; j++ ) {
	  map[ j + n ].p[k] = k ;
	}
      } else {
	if( sim_params[k] == true ) {
	  for( j = 0 ; j < Ndata[i] ; j++ ) {
	    map[ j + n ].p[k] = k ;
	  }
	} else {
	  for( j = 0 ; j < Ndata[i] ; j++ ) {
	    map[ j + n ].p[k] = addition ;
	  }
	  addition++ ;
	}
      }
    }
    n += Ndata[i] ;
    //
  }

#ifdef verbose
  // loop this stuff
  n = 0 ;
  for( i = 0 ; i < Nsims ; i++ ) {
    for( j = 0 ; j < Ndata[i] ; j++ ) {
      printf( "THIS -> %zu :: " , n ) ;
      //
      for( k = 0 ; k < Nparams ; k++ ) {
	printf( " %zu " , map[ n ].p[ k ] ) ; 
      }
      printf( "\n" ) ;
      //
      n++ ;
    }
  }
#endif
  
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
