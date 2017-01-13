/**
   @file poly_coefficients.c
   @brief compute polynomial coefficients using SVD

   Using GSL's SVD we attempt to compute the coefficients of a polynomial

   i.e. solve the (possibly) over-constrained VanderMonde equation


   | 1 x_0 x_0^2 x_0^3 .... x_0^N | | a_0 |   | y( x_0 ) |
   | 1 x_1 x_1^2 x_1^3 .... x_1^N | | a_1 | = | y( x_1 ) |
   | ............................ | | ... | = |  .....   |
   | 1 .................... x_M^N | | a_N |   | y( x_M ) |
 
   Where M >= N and the a's are the coefficients of the polynomial
 */
#include "gens.h"
#include "svd.h"

int
compute_coefficients( double *coeffs ,
		      double *chisq ,
		      const double *y ,
		      const double *sigma ,
		      const double *x ,
		      const size_t M ,
		      const size_t N )
{
  if( M == 0 || N == 0 ) return FAILURE ;
  
  // allocate the matrices
  double **A = malloc( M * sizeof( double* ) ) ;
  double **Ainv = malloc( N * sizeof( double* ) ) ;
  register double sum , x0 ;
  size_t i , j , FLAG = 0 ;

  // matrix allocations
  for( i = 0 ; i < M ; i++ ) {
    A[i] = (double*)malloc( N * sizeof( double ) ) ;
    if( i < N ) 
      Ainv[i] = (double*)malloc( M * sizeof( double ) ) ;
  }

  // allocate the matrix with increasing powers of x
  for( i = 0 ; i < M ; i++ ) {
    x0 = 1 ;
    for( j = 0 ; j < N ; j++ ) {
      A[i][j] = x0 ;
      x0 *= x[i] ;
    }
  }

  // if the SVD screws up we set failure and free allocations
  if( svd_inverse( Ainv , (const double**)A , N , M , 1E-8 , true
		   ) == FAILURE ) {
    FLAG = FAILURE ;
    goto FREE ;
  }

  // multiply the inverse by the data to obtain the coefficients
  for( i = 0 ; i < N ; i++ ) {
    sum = 0.0 ;
    for( j = 0 ; j < M ; j++ ) {
      sum += Ainv[i][j] * y[j] ;
    }
    coeffs[i] = sum ;
  }

  // compute residual r += ( y[ i ] - A * c )
  *chisq = 0.0 ;
  register double ri ;
  for( j = 0 ; j < M ; j++ ) {
    register double res = 0.0 ;
    for( i = 0 ; i < N ; i++ ) {
      res += A[j][i] * coeffs[i] ;
    }
    ri = ( y[ j ] - res ) / ( sigma[ j ] ) ;
    *chisq += ri * ri ;
  }

 FREE :
  // free this stuff
  for( i = 0 ; i < M ; i++ ) {
    free( A[i] ) ;
    if( i < N ) {
      free( Ainv[i] ) ;
    }
  }
  free( A ) ;
  free( Ainv ) ;
  return FLAG ;
}

// just a helpful little routine
void
write_polynomial( const double *__restrict coeffs ,
		  const size_t POLY_ORDER )
{
  size_t j ;
  printf( "y = " ) ;
  for( j = 0 ; j < POLY_ORDER-1 ; j++ ) {
    printf( "%f * x^{%zu} + " , coeffs[j] , j ) ;
  }
  printf( "%f * x^{%zu) \n" , coeffs[j] , j ) ;
  printf( "\n" ) ;
}
