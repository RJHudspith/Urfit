/**
   @file Ren_rats.c
   @brief renormalise ratios and convert to SUSY basis
 */
#include "gens.h"

#include "resampled_ops.h"

int
renormalise_rats( struct input_params *Input )
{
  // change to 7 when I can be bothered so we can multiply by m_K^2/f_K^2
  if( Input -> Data.Nsim != 5 ) {
    fprintf( stderr , "Code is very specific about number of inputs\n" ) ;
    return FAILURE ;
  }

  if( Input -> Data.Ndata[0] != 25 ) {
    fprintf( stderr , "Zmatrix should be first and have 25 elements\n" ) ;
    return FAILURE ;
  }

  // do the renormalisation
  struct resampled *Q = malloc( 5 * sizeof( struct resampled ) ) ;
  struct resampled *O = malloc( 5 * sizeof( struct resampled ) ) ;
  size_t i , j ;
  for( i = 0 ; i < 5 ; i++ ) {
    Q[i] = init_dist( NULL , Input -> Data.y[25].NSAMPLES ,
		      Input -> Data.y[25].restype ) ;
    O[i] = init_dist( NULL , Input -> Data.y[25].NSAMPLES ,
		      Input -> Data.y[25].restype ) ;
  }

  struct resampled temp = init_dist( NULL , Input -> Data.y[25].NSAMPLES ,
				     Input -> Data.y[25].restype ) ;
  for( i = 1 ; i < 5 ; i++ ) {
    for( j = 1 ; j < 5 ; j++ ) {
      equate( &temp , Input -> Data.y[25+j-1] ) ;
      mult( &temp , Input -> Data.y[j+i*5] ) ;
      add( &Q[i] , temp ) ; 
    }
    divide( &Q[i] , Input -> Data.y[0] ) ;

    printf( "REN %f %f \n" , Q[i].avg , Q[i].err ) ;
  }

  // convert to SUSY basis
  equate( &O[1] , Q[3] ) ;
  
  equate( &temp , Q[3] ) ;
  subtract( &temp , Q[4] ) ;
  mult_constant( &temp , -0.5 ) ;
  equate( &O[2] , temp ) ;

  equate( &O[3] , Q[2] ) ;

  equate( &temp , Q[1] ) ;
  mult_constant( &temp , -0.5 ) ;
  equate( &O[4] , temp ) ;

  for( i = 1 ; i < 5 ; i++ ) {
    printf( "SUSY_%zu %f %f \n" , i , O[i].avg , O[i].err ) ;
  }

  for( i = 0 ; i < 5 ; i++ ) {
    free( O[i].resampled ) ;
    free( Q[i].resampled ) ;
  }
  free( O ) ;
  free( Q ) ;

  free( temp.resampled ) ;

  return SUCCESS ;
}
