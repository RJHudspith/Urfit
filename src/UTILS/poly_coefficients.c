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
  if( svd_inverse( Ainv , (const double**)A , N , M , 1E-14 , true
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

// computes using ( U U^T ) U y
int
compute_coefficients2( double *coeffs ,
		      double *chisq ,
		      const double *y ,
		      const double *sigma ,
		      const double *x ,
		      const size_t M ,
		      const size_t N )
{
  if( M == 0 || N == 0 ) return FAILURE ;

  // solve by LU? rolls back to column-balanced SVD if it can't
  gsl_matrix *alpha = gsl_matrix_alloc( N , N ) ;
  gsl_vector *beta  = gsl_vector_alloc( N ) ;
  gsl_vector *delta = gsl_vector_alloc( N ) ;
  gsl_permutation *perm = gsl_permutation_alloc( N ) ;
  
  register double sum , x0 ;
  size_t i , j , k , FLAG = 0 ;

  // matrix allocations
  double **A = malloc( M * sizeof( double* ) ) ;
  for( i = 0 ; i < M ; i++ ) {
    A[i] = malloc( N * sizeof( double ) ) ;
  }

  // allocate the matrix with increasing powers of x
  for( i = 0 ; i < M ; i++ ) {
    x0 = 1 ;
    for( j = 0 ; j < N ; j++ ) {
      A[i][j] = x0 ;
      x0 *= x[i] ;
    }
  }

  // initialise alpha as the square of A 
  for( i = 0 ; i < N ; i++ ) {
    for( j = i ; j < N ; j++ ) {
      register double sum = 0.0 ;
      for( k = 0 ; k < M ; k++ ) {
	sum += A[k][i] * A[k][j] ;
      }
      gsl_matrix_set( alpha , i , j , sum ) ;
      gsl_matrix_set( alpha , j , i , sum ) ;
    }
  }

  printf( "\n" ) ;
  for( i = 0 ; i < N ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      printf( " %1.2f " , gsl_matrix_get( alpha , i , j ) ) ;
    }
    printf( "\n" ) ;
  }

  // compute beta
  for( i = 0 ; i < N ; i++ ) {
    register double sum = 0.0 ;
    for( j = 0 ; j < M ; j++ ) {
      sum += A[j][i] * y[j] ;
    }
    gsl_vector_set( beta , i , sum ) ;
  }
  
  // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
  bool svd_flag = false ;
  int signum ;
  if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
    fprintf( stderr , "[POLY_COEFF] LU decomp failed, rolling back to SVD\n" ) ;
    svd_flag = true ;
  }
  if( svd_flag != true ) {
    if( gsl_linalg_LU_solve( alpha , perm , beta , delta ) != GSL_SUCCESS ) {
      fprintf( stderr , "[POLY_COEFF] LU solve failure \n" ) ;
      goto FREE ;
    }
    // set the coefficients
    for( i = 0 ; i < N ; i++ ) {
      coeffs[i] = gsl_vector_get( delta , i ) ;
    }
  } else {
    printf( "[POLY] regressing to the SVD \n" ) ;
    
    double **Ainv = malloc( N * sizeof( double* ) ) ;
    double **Anew = malloc( N * sizeof( double* ) ) ;
    for( i = 0 ; i < N ; i++ ) {
      Ainv[i] = malloc( N * sizeof( double ) ) ;
      Anew[i] = malloc( N * sizeof( double ) ) ;
      for( j = 0 ; j < N ; j++ ) {
	Anew[i][j] = gsl_matrix_get( alpha , i , j ) ;
      }
    }
    // if the SVD screws up we set failure and free allocations
    if( svd_inverse( Ainv , (const double**)Anew , N , N , 1E-14 , true
		     ) == FAILURE ) {
      FLAG = FAILURE ;
      goto FREE ;
    }
    
    // multiply the inverse by the data to obtain the coefficients
    for( i = 0 ; i < N ; i++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < N ; j++ ) {
	sum += Ainv[i][j] * gsl_vector_get( beta , j ) ;
      }
      coeffs[i] = sum ;
    }

    // memfree
    for( i = 0 ; i < M ; i++ ) {
      free( Ainv[i] ) ;
      free( Anew[i] ) ;
    }
    free( Ainv ) ;
    free( Anew ) ;
  }

  // compute residual r += ( y[ i ] - A * c )
  *chisq = 0.0 ;
  register double ri ;
  for( j = 0 ; j < M ; j++ ) {
    register double res = 0.0 ;
    for( i = 0 ; i < N ; i++ ) {
      res += A[j][i] * coeffs[i] ;
    }
    ri = ( y[ j ] - res ) ;
    *chisq += ri * ri ;
  }

 FREE :

  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  // free this stuff
  for( i = 0 ; i < M ; i++ ) {
    free( A[i] ) ;
  }
  free( A ) ;
  
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
