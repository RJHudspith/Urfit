/**
   @file pade_coefficients.c
   @brief given polynomial coefficients, compute the pade representation


   p( x ) = \pi_0 + ( \pi_1 x + ... + \pi_n x^n ) / ( 1 + ... + \pi_{n+m}x^{m} )
 */
#include <stdlib.h>
#include "svd.h"

#define SUCCESS (!FAILURE)
#define FAILURE (0)

// pade coefficient function
int
pades( double *pade_coeffs ,
       const double *poly_coeffs ,
       const size_t n ,
       const size_t m )
{
  // is guaranteed
  pade_coeffs[ 0 ] = poly_coeffs[ 0 ] ;

  // I do the simple cases here 
  switch( n ) {
  case 1 :
    switch( m ) {
    case 1 :     // ( 1 , 1 ) pade
      pade_coeffs[ 1 ] =  poly_coeffs[ 1 ] ;
      pade_coeffs[ 2 ] = -poly_coeffs[ 2 ] / poly_coeffs[ 1 ] ;
      return SUCCESS ;
    case 2 :     // ( 1 , 2 ) pade
      pade_coeffs[ 1 ] =  poly_coeffs[ 1 ] ;
      pade_coeffs[ 2 ] = -poly_coeffs[ 2 ] / poly_coeffs[ 1 ] ;
      pade_coeffs[ 3 ] = ( poly_coeffs[ 2 ] * poly_coeffs[ 2 ] - 
			   poly_coeffs[ 1 ] * poly_coeffs[ 3 ] ) / 
	( poly_coeffs[1] * poly_coeffs[1] ) ;
      return SUCCESS ;
    }
    break ;
  case 2 :
    switch( m ) {
    case 1 :     // ( 2 , 1 ) pade
      pade_coeffs[ 1 ] = poly_coeffs[ 1 ] ;
      pade_coeffs[ 2 ] = poly_coeffs[ 2 ] - 
	( poly_coeffs[ 3 ] * poly_coeffs[ 1 ] ) / poly_coeffs[ 2 ] ;
      pade_coeffs[ 3 ] = -poly_coeffs[ 3 ] / poly_coeffs[ 2 ] ;
      return SUCCESS ;
    }
    break ;
  }

  // the matrix order equation we must solve
  const size_t order = n + m + 1 ;

  // initialise solve matrix and its inverse
  double **sub_solve = malloc( m * sizeof( double* ) ) ;
  double **Ainv = malloc( m * sizeof( double* ) ) ;

  // vectors we use
  double *solves = malloc( m * sizeof( double ) ) ;
  double *left_side = malloc( ( order - 1 ) * sizeof( double ) ) ;
  double *right_side = malloc( ( order - 1 ) * sizeof( double ) ) ;

  // set up the submatrix we need to solve for the
  // left hand side
  size_t i , j , FLAG = FAILURE ;
  for( i = 0 ; i < m ; i++ ) {
    sub_solve[ i ] = malloc( m * sizeof( double) ) ;
    Ainv[ i ] = malloc( m * sizeof( double) ) ;
    for( j = 0 ; j < m ; j++ ) {
      if( i - j + n > 0 ) {
	sub_solve[i][j] = poly_coeffs[ i - j + n ] ;
      } else {
	sub_solve[i][j] = 0.0 ;
      }
    }
    solves[i] = -poly_coeffs[ i + n + 1 ] ;
    #ifdef VERBOSE
    printf( "sols [ %d ] = %d :: %f \n" , i , m + i + 1 , solves[i] ) ;
    #endif
  }

  // printhimout
#ifdef VERBOSE
  printmatrix( (const double**)sub_solve , m , m ) ;
#endif

  // perform an inverse
  if( svd_inverse( Ainv , (const double**)sub_solve , m , m ) == FAILURE ) {
    FLAG = FAILURE ;
    goto FREE ;
  }

  // multiply ainverse by sols to get left hand side coefficients
  left_side[ 0 ] = 1.0 ;
  for( i = 1 ; i < order - 1 ; i++ ) {
    left_side[ i ] = 0.0 ;
    if( i <= m ) {
      for( j = 0 ; j < m ; j++ ) {
	left_side[ i ] += Ainv[i-1][j] * solves[j] ;
      }
    }
  }

  // multiply coefficient matrix by left hand side
  for( i = 0 ; i < order-1 ; i++ ) {
    register double sum = poly_coeffs[ i + 1 ] ;
    for( j = 1 ; j <= i ; j++ ) {
      sum += poly_coeffs[ i - j + 1 ] * left_side[ j ] ;
    }
    right_side[ i ] = sum ;
    #ifdef VERBOSE
    printf( "\n" ) ;
    printf( "right[ %d ] = %f \n" , i , right_side[ i ] ) ;
    #endif
  }

  // numerator is first in our scheme
  for( i = 0 ; i < n-1 ; i++ ) {
    pade_coeffs[ i + 1 ] = right_side[ i ] ;
  }
  pade_coeffs[ i + 1 ] = right_side[ i ] ;
  // denominator is last
  for( i = 0 ; i < m-1 ; i++ ) {
    pade_coeffs[ i + n + 1 ] = left_side[ i + 1 ] ;
  } 
  pade_coeffs[ i + n + 1 ] = left_side[ i + 1 ] ;

  // free workspaces
 FREE :
  for( i = 0 ; i < m ; i++ ) {
    free( sub_solve[i] ) ;
    free( Ainv[i] ) ;
  }
  free( sub_solve ) ;
  free( Ainv ) ;
  free( solves ) ;
  free( left_side ) ;
  free( right_side ) ;

  return FLAG ;
}
