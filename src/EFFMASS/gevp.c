/**
   @file gevp.c
   @brief solves a generalised eigenvalue problem

   Solves the system

   A.v = \lambda B.v
   
   where A and B are real matrices
 */
#include <stdio.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define FAILURE -1
#define SUCCESS !FAILURE

// uses GSL to solve a generalised eigenvalue problem
// returns the real part of the eigenvalues
// solves A.v = \lambda B.v
// where A and B are real matrices
int
solve_gevp( double *re_evalues , 
	    const double *A , // linearised nxn matrix
	    const double *B , // linearised nxn matrix
	    const size_t n ) 
{
  int flag = SUCCESS ;

  // allocations for GSL
  gsl_eigen_gen_workspace *work = gsl_eigen_gen_alloc ( n ) ;
  gsl_matrix *a  = gsl_matrix_alloc( n , n ) ;
  gsl_matrix *b  = gsl_matrix_alloc( n , n ) ;
  gsl_vector_complex *alpha  = gsl_vector_complex_alloc( n ) ;
  gsl_vector *beta  = gsl_vector_alloc( n ) ;

  // set the matrices
  size_t i , j ;
  for( i = 0 ; i < n ; i++ ){
    for( j = 0 ; j < n ; j++ ) {
      gsl_matrix_set( a , i , j , A[ j + i*n ] ) ;
      gsl_matrix_set( b , i , j , B[ j + i*n ] ) ;
    }
  }

  // perform decomposition
  const int err = gsl_eigen_gen( a , b , alpha , beta , work ) ;
  if( err != 0 ) {
    printf( "%s\n" , gsl_strerror( err ) ) ;
    printf( "Aborting\n" ) ;
    flag = FAILURE ;
    goto free ;
  }

  // get the real and imaginary parts
  for( i = 0 ; i < n ; i++ ) {
    re_evalues[ i ] = ( gsl_vector_complex_get( alpha , i ).dat[0] ) 
      / gsl_vector_get( beta , i ) ;
  }

  // memfreeze
 free :
  gsl_matrix_free( a ) ;
  gsl_matrix_free( b ) ;
  gsl_vector_complex_free( alpha ) ;
  gsl_vector_free( beta ) ;
  gsl_eigen_gen_free( work ) ;

  return flag ;
}
