/**
   @file gevp.c
   @brief solves a generalised eigenvalue problem

   Solves the system

   A.v = \lambda B.v
   
   where A and B are real matrices
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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
	    const size_t n ,
	    const bool before ,
	    const bool write_evalues ) 
{
  int flag = SUCCESS ;

  // allocations for GSL
  gsl_eigen_genv_workspace *work = gsl_eigen_genv_alloc( n ) ;

  gsl_matrix *a  = gsl_matrix_alloc( n , n ) ;
  gsl_matrix *b  = gsl_matrix_alloc( n , n ) ;

  gsl_vector_complex *alpha  = gsl_vector_complex_alloc( n ) ;
  gsl_vector *beta  = gsl_vector_alloc( n ) ;
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc( n , n ) ;

  // set the matrices
  size_t i , j ;
  for( i = 0 ; i < n ; i++ ){
    for( j = 0 ; j < n ; j++ ) {
      gsl_matrix_set( a , i , j , A[ j + i*n ] ) ;
      gsl_matrix_set( b , i , j , B[ j + i*n ] ) ;
    }
  }

  // perform decomposition
  const int err = gsl_eigen_genv( a , b , alpha , beta , evec , work ) ;

  if( err != 0 ) {
    printf( "%s\n" , gsl_strerror( err ) ) ;
    printf( "Aborting\n" ) ;
    flag = FAILURE ;
    goto free ;
  }

  if( before ) {
    gsl_eigen_genv_sort( alpha , beta , evec , GSL_EIGEN_SORT_ABS_ASC ) ;
  } else {
    gsl_eigen_genv_sort( alpha , beta , evec , GSL_EIGEN_SORT_ABS_DESC ) ;
  }
  
  // get the real and imaginary parts
  for( i = 0 ; i < n ; i++ ) {
    // if we want them written
    if( write_evalues == true ) {
      fprintf( stdout , "EVALUE_%zu %e %e \n" ,
	       i , ( gsl_vector_complex_get( alpha , i ).dat[0] ) / gsl_vector_get( beta , i ) ,
	       ( gsl_vector_complex_get( alpha , i ).dat[1] ) / gsl_vector_get( beta , i ) ) ;
    }

    re_evalues[ i ] = ( gsl_vector_complex_get( alpha , i ).dat[0] ) 
      / gsl_vector_get( beta , i ) ;
  }
  
  // memfreeze
 free :
  gsl_matrix_free( a ) ;
  gsl_matrix_free( b ) ;

  gsl_matrix_complex_free( evec ) ;
  
  gsl_vector_complex_free( alpha ) ;
  gsl_vector_free( beta ) ;
  gsl_eigen_genv_free( work ) ;

  return flag ;
}
