/**
   @file inverse.c
   @brief given some distributions it inverts them
 */
#include "gens.h"

//#define BD

int
fit_inverse( struct input_params *Input )
{
  size_t i , j , shift = 0 ;

  if( Input -> Data.Nsim != 25 ) {
    printf( "Nsim not as expected %zu \n" , Input -> Data.Nsim ) ;
    return FAILURE ;
  }
  
  gsl_matrix *alpha = gsl_matrix_alloc( 5 , 5 ) ;
  gsl_matrix *ainv  = gsl_matrix_alloc( 5 , 5 ) ;
  gsl_permutation *perm = gsl_permutation_alloc( 5 ) ;
  int signum ;
  
  // invert each case by LU
  size_t l , k ;
  printf( "Zinverse matrix\n" ) ;
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {
      printf( "{%e %e} " , Input -> Data.y[ k + 5*l ].avg ,
	      Input -> Data.y[ k + 5*l ].err ) ;
    }
    printf( "\n" ) ;
  }
  printf( "\n" ) ;
  
  for( i = 0 ; i < Input -> Data.y[0].NSAMPLES ; i++ ) {

    #ifdef BD
    for( l = 0 ; l < 5 ; l++ ) {
      for( k = 0 ; k < 5 ; k++ ) {
	gsl_matrix_set( alpha , l , k , 0.0 ) ;
      }
    }
    gsl_matrix_set( alpha , 0 , 0 , Input -> Data.y[0].resampled[i] ) ;
    
    gsl_matrix_set( alpha , 1 , 1 , Input -> Data.y[6].resampled[i] ) ;
    gsl_matrix_set( alpha , 1 , 2 , Input -> Data.y[7].resampled[i] ) ;
    gsl_matrix_set( alpha , 2 , 1 , Input -> Data.y[11].resampled[i] ) ;
    gsl_matrix_set( alpha , 2 , 2 , Input -> Data.y[12].resampled[i] ) ;

    gsl_matrix_set( alpha , 3 , 3 , Input -> Data.y[18].resampled[i] ) ;
    gsl_matrix_set( alpha , 3 , 4 , Input -> Data.y[19].resampled[i] ) ;
    gsl_matrix_set( alpha , 4 , 3 , Input -> Data.y[23].resampled[i] ) ;
    gsl_matrix_set( alpha , 4 , 4 , Input -> Data.y[24].resampled[i] ) ;
    
    #else
    for( l = 0 ; l < 5 ; l++ ) {
      for( k = 0 ; k < 5 ; k++ ) {
	gsl_matrix_set( alpha , l , k ,
			Input -> Data.y[ k + 5*l ].resampled[i] ) ;
      }
    }
    #endif
    // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
    if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
      printf( "[LM] LU decomp broken?\n" ) ;
      return !GSL_SUCCESS ;
    }
    if( gsl_linalg_LU_invert( alpha , perm , ainv ) != GSL_SUCCESS ) {
      printf( "[LM] LU solve broken?\n" ) ;
      return !GSL_SUCCESS ;
    }
    // poke it into the x direction
    for( l = 0 ; l < 5 ; l++ ) {
      for( k = 0 ; k < 5 ; k++ ) {
	Input -> Data.x[ k + 5*l ].resampled[i] =
	  gsl_matrix_get( ainv , l , k ) ;
      }
    }
    //
  }

  printf( "\nZ matrix\n" ) ;
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {
      compute_err( &Input -> Data.x[ k + 5*l ] ) ;
      printf( "%f(%f) & " , l , k ,
	      Input -> Data.x[ k + 5*l ].avg ,
	      Input -> Data.x[ k + 5*l ].err ) ;
    }
    printf( "\\\\ \n" ) ;
  }
				   
  return SUCCESS ;
}
