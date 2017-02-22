/**
   @file GLS_pade.c
   @brief generalised least squares for pade functions

   Solves Z \beta = \alpha for \beta and \alpha

   representation is

   y = A \alpha / A \beta

   where A is a Vandermonde matrix and \beta[0] = 1
 */
#include "gens.h"

#include "chisq.h"
#include "svd.h"

#ifdef VERBOSE
static void
printmatrix( const double **mat ,
	     const int NROWS , 
	     const int NCOLS )
{
  printf( "\n" ) ;
  int i , j ;
  for( i = 0 ; i < NROWS ; i++ ) {
    printf( "|" ) ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      printf( " %e " , mat[i][j] ) ;
    }
    printf( "|\n" ) ;
  }
  printf( "\n" ) ;
  return ;
}
#endif

// pade coefficient function
int
pades( double *pade_coeffs ,
       const double **Z ,
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
    for( j = 0 ; j < ( m ) ; j++ ) {
      gsl_matrix_set( alpha , i , j , Z[ i + n ][ j + 1 ] ) ;
      printf( "SUB %zu %zu -> %e\n" , i , j , Z[i+n][j+1] ) ;
    }
    gsl_vector_set( beta , i , -Z[ n + i ][ 0 ] ) ;
  }

  // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
  int signum , Flag = SUCCESS ;
  if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
    fprintf( stderr , "[POLY_COEFF] LU decomp failed, rolling back to SVD\n" ) ;
    Flag = FAILURE ;
  }
  if( gsl_linalg_LU_solve( alpha , perm , beta , delta ) != GSL_SUCCESS ) {
    fprintf( stderr , "[POLY_COEFF] LU solve failure \n" ) ;
    Flag = FAILURE ;
  }

  // multiply ainverse by sols to get left hand side coefficients
  for( j = 0 ; j < m ; j++ ) {
    pade_coeffs[ n + i - 1 ] = gsl_vector_get( delta , j ) ;
  }

  // multiply coefficient matrix by left hand side
  for( i = 0 ; i < n ; i++ ) {
    register double sum = Z[ i ][ 0 ] ;
    for( j = 0 ; j < m ; j++ ) {
      sum += Z[ i ][ j + 1 ] * pade_coeffs[ n + j ] ;
    }
    pade_coeffs[ i ] = sum ;
    #ifdef VERBOSE
    printf( "\n" ) ;
    printf( "right[ %d ] = %f \n" , i , pade_coeffs[ i ] ) ;
    #endif
  }

#ifdef VERBOSE
  // numerator is first in our scheme
  printf( "( (%1.2f) q^{%d} " , pade_coeffs[ 0 ] , 2*(0) ) ;
  for( i = 1 ; i < n ; i++ ) {
    printf( " + (%1.2f) q^{%d} " , pade_coeffs[ i ] , 2*(i) ) ;
  }
  printf( ") / ( 1 " ) ;

  // denominator is last
  for( i = 0 ; i < m ; i++ ) {
    printf( "+ (%1.2f) q^{%d} " , pade_coeffs[ i + n ] , 2*(i+1) ) ;
  }
  printf( " ) \n" ) ;
#endif

  // free the gsl vectors
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  return FLAG ;
}

// computes the Z matrix
static double **
compute_newleft( const double *xdata ,
		 const double *ydata ,
		 const double **W ,
		 const corrtype CORRFIT ,
		 const size_t NROWS , // this is the number of data
		 const size_t NCOLS ) // this is the order of the pade
{
  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( NCOLS , NCOLS ) ;
  gsl_matrix *alinv = gsl_matrix_alloc( NCOLS , NCOLS ) ;
  gsl_permutation *perm = gsl_permutation_alloc( NCOLS ) ;
  int signum ;
  
  double **A = malloc( NROWS * sizeof( double* ) ) ;

  // compute the matrix A is the vandermonde matrix
  size_t i , j , k , l ;
  for( i = 0 ; i < NROWS ; i++ ) {
    A[i] = malloc( NCOLS * sizeof( double ) ) ;
    register double xs = 1 , xx = xdata[i] ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      A[ i ][ j ] = xs ;
      xs *= ( xx ) ;
    }
  }

#ifdef VERBOSE
  printf( "A MATRIX \n" ) ;
  printmatrix( A , NROWS , NCOLS ) ;
#endif
  
  // Z == A^T * Cov * A
  for( i = 0 ; i < NCOLS ; i++ ) {
    for( j = 0 ; j < NCOLS ; j++ ) {
      register double sum = 0.0 ;
      for( k = 0 ; k < NROWS ; k++ ) {
	switch( CORRFIT ) {
	case UNWEIGHTED :
	  sum += A[k][i] * A[k][j] ;
	  break ;
	case UNCORRELATED :
	  sum += A[k][i] * W[0][k] * A[k][j] ;
	  break ;
	case CORRELATED :
	  for( l = 0 ; l < NROWS ; l++ ) {
	    sum += A[k][i] * W[k][l] * A[l][j] ;
	  }
	  break ;
	}
      }
      gsl_matrix_set( alpha , i , j , sum ) ;
    }
  }

  if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
    fprintf( stderr , "[PADE_COEFF] LU decomp failed, rolling back to SVD\n" ) ;
    return NULL ;
  }
  
  if( gsl_linalg_LU_invert( alpha , perm , alinv ) != GSL_SUCCESS ) {
    fprintf( stderr , "[PADE_COEFF] LU invert failure \n" ) ;
    return NULL ;
  }

  // compute the right hand part y' = ( A^T Cov y^T A )
  for( i = 0 ; i < NCOLS ; i++ ) {
    for( j = 0 ; j < NCOLS ; j++ ) {
      register double sum = 0.0 ;
      for( k = 0 ; k < NROWS ; k++ ) {
	switch( CORRFIT ) {
	case UNWEIGHTED :
	  sum += A[k][i] * ydata[k] * A[k][j] ;
	  break ;
	case UNCORRELATED :
	  sum += A[k][i] * W[0][k] * ydata[k] * A[k][j] ;
	  break ;
	case CORRELATED :
	  for( l = 0 ; l < NROWS ; l++ ) {
	    sum += A[k][i] * W[k][l] * ydata[l] * A[l][j] ;
	  }
	  break ;
	}
      }
      gsl_matrix_set( alpha , i , j , sum ) ;
    }
  }

#ifdef VERBOSE
  printf( "Z2 set \n" ) ;
  printmatrix( Z , NCOLS , NCOLS ) ;
#endif

  // reallocate A
  A = realloc( A , NCOLS * sizeof( double * ) ) ;

  // finally multiply the two matrices
  for( i = 0 ; i < NCOLS ; i++ ) {
    for( j = 0 ; j < NCOLS ; j++ ) {
      register double sum = 0.0 ;
      for( k = 0 ; k < NCOLS ; k++ ) {
	sum += gsl_matrix_get( alinv , i , k ) *
	  gsl_matrix_get( alpha , k , j ) ;
      }
      A[i][j] = sum ;
    }
  }

  #ifdef VERBOSE
  printmatrix( A , NCOLS , NCOLS ) ;
  #endif

  // free the gsl vectors
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( alinv ) ;
  gsl_permutation_free( perm ) ;

  return A ;
}

// perform a generalised least squares pade fit iteration
int
gls_pade_iter( void *fdesc ,
	       const void *data ,
	       const double **W ,
	       const double TOL )
{
  // point to the fit function
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;

    // point at the data
  struct data *Data = (struct data*)data ;

  // set the priors
  Fit -> f.Prior = Fit -> Prior ;

  // number of poly coefficients is Nlogic
  const size_t N = Fit -> f.N ;
  const size_t M = Fit -> Nlogic ;

  double **Z = compute_newleft( Data -> x , Data -> y , W ,
				Fit -> f.CORRFIT , N , M ) ;

  pades( Fit -> f.fparams , (const double**)Z , Data -> N , Data -> M ) ;

  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

  size_t i ;
  for( i = 0 ; i < M ; i++ ) {
    free( Z[i] ) ;
  }
  free( Z ) ;

  return SUCCESS ;
}
