/**
   @file pade_coefficients.c
   @brief given polynomial coefficients, compute the pade representation


   p( x ) = ( \pi_0 + \pi_1 x + ... + \pi_n x^n ) / ( 1 + ... + \pi_{n+m}x^{m} )
 */
#include "gens.h"

// compute pade coefficients from a polynomial series
int
pades_from_poly( double *pade_coeffs ,
		 const double *poly_coeffs ,
		 const size_t n ,
		 const size_t m )
{
  // the matrix order equation we must solve
  size_t i , j ;
  for( i = 0 ; i < n + m ; i++ ) {
    pade_coeffs[i] = 0.0 ;
  }

  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( m , m ) ;
  gsl_vector *beta  = gsl_vector_alloc( m ) ;
  gsl_vector *delta = gsl_vector_alloc( m ) ;
  gsl_permutation *perm = gsl_permutation_alloc( m ) ;

  // flag for whether it worked
  int FLAG = FAILURE ;

  // set up the submatrix we need to solve for the
  // left hand side
  for( i = 0 ; i < m ; i++ ) {
    for( j = 0 ; j < m ; j++ ) {
      if( ( i + n ) < ( j + 1 ) ) {
	gsl_matrix_set( alpha , i , j , 0 ) ;
      } else {
	gsl_matrix_set( alpha , i , j , poly_coeffs[ i + n - j - 1 ] ) ;
      }
    }
    gsl_vector_set( beta , i , -poly_coeffs[ i + n ] ) ;
  }

  // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
  int signum , Flag = SUCCESS ;
  if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
    fprintf( stderr , "[POLY_COEFF] LU decomp failed \n" ) ;
    Flag = FAILURE ;
  }
  if( gsl_linalg_LU_solve( alpha , perm , beta , delta ) != GSL_SUCCESS ) {
    fprintf( stderr , "[POLY_COEFF] LU solve failure \n" ) ;
    Flag = FAILURE ;
  }

  // multiply ainverse by sols to get left hand side coefficients
  for( j = 0 ; j < m ; j++ ) {
    pade_coeffs[ n + j ] = gsl_vector_get( delta , j ) ;
  }

  // multiply coefficient matrix by left hand side
  for( i = 0 ; i < n ; i++ ) {
    register double sum = poly_coeffs[ i ] ;
    for( j = 1 ; j <= i ; j++ ) {
      #ifdef VERBOSE
      fprintf( stdout , "Coeffs :: (%zu,%zu) %zu -> %zu \n" ,
	       i , j , i-j , n+j-1 ) ;
      #endif
      sum += poly_coeffs[ i - j  ] * pade_coeffs[ n + j - 1 ] ;
    }
    pade_coeffs[ i ] = sum ;
    #ifdef VERBOSE
    fprintf( stdout , "\n" ) ;
    fprintf( stdout , "right[ %zu ] = %f \n" , i , pade_coeffs[ i ] ) ;
    #endif
  }

  #ifdef VERBOSE
  // numerator is first in our scheme
  fprintf( stdout , "( (%1.2f) q^{%d} " , pade_coeffs[ 0 ] , 2*(0) ) ;
  for( i = 1 ; i < n ; i++ ) {
    fprintf( stdout , " + (%1.2f) q^{%zu} " , pade_coeffs[ i ] , 2*(i) ) ;
  }
  fprintf( stdout , ") / ( 1 " ) ;

  // denominator is last
  for( i = 0 ; i < m ; i++ ) {
    fprintf( stdout , "+ (%1.2f) q^{%zu} " , pade_coeffs[ i + n ] , 2*(i+1) ) ;
  }
  fprintf( stdout , " ) \n" ) ;
  #endif

  // free the gsl vectors
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;
  
  return FLAG ;
}
