/**
   @file gevp.c
   @brief solves a generalised eigenvalue problem

   Solves the system

   A.v = \lambda B.v
   
   where A and B are real matrices
 */
#include "gens.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
/*
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
#define FAILURE -1
#define SUCCESS !FAILURE
*/

// uses GSL to solve a generalised eigenvalue problem
// returns the real part of the eigenvalues
// solves A.v = \lambda B.v
// where A and B are real matrices
static int
gevp( gsl_matrix *a ,
      gsl_matrix *b ,
      gsl_eigen_genv_workspace *work ,
      gsl_vector_complex *alpha ,
      gsl_vector *beta ,
      gsl_matrix_complex *evec ,
      const size_t n ,
      const bool before ,
      const bool write_evalues ) 
{
  size_t i , j ;
  int flag = SUCCESS ;
  
  // perform decomposition
  const int err = gsl_eigen_genv( a , b , alpha , beta , evec , work ) ;

  if( err != 0 ) {
    printf( "%s\n" , gsl_strerror( err ) ) ;
    printf( "Aborting\n" ) ;
    flag = FAILURE ;
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
	       i , ( gsl_vector_complex_get( alpha , i ).dat[0] )
	       / gsl_vector_get( beta , i ) ,
	       ( gsl_vector_complex_get( alpha , i ).dat[1] )
	       / gsl_vector_get( beta , i ) ) ;
    }
  }

  return flag ;
}

// take the SVD
static int 
svd_TLS( gsl_matrix *a ,
	 gsl_matrix *b ,
	 const double *C0 ,
	 const double *C1 ,
	 const size_t NROWS ,
	 const size_t NCOLS )
{
  // allocations
  gsl_matrix *U  = gsl_matrix_alloc( 2*NROWS , NCOLS ) ;
  gsl_matrix *Q  = gsl_matrix_alloc( NCOLS , NCOLS ) ;
  gsl_vector *S  = gsl_vector_alloc( NCOLS ) ;
  gsl_vector *WORK  = gsl_vector_alloc( NCOLS ) ;

  size_t i , j ;
  int FLAG = SUCCESS ;

  // initialize our matrix
  for( i = 0 ; i < NROWS ; i++ ){
    for( j = 0 ; j < NCOLS ; j++ ) {
      gsl_matrix_set( U , i , j , C0[j+i*NCOLS] ) ;
      gsl_matrix_set( U , i + NROWS , j , C1[j+i*NCOLS] ) ;
    }
  }

  // perform the decomposition
  if( gsl_linalg_SV_decomp_jacobi( U , Q , S ) ) {
    printf( "[SVD] GSL SVD comp failure \n" ) ;
    FLAG = FAILURE ;
    goto FREE ;
  }

  for( i = 0 ; i < NCOLS ; i++ ) {
    for( j = 0 ; j < NCOLS ; j++ ) {
      gsl_matrix_set( a , i , j , gsl_matrix_get( U , i , j ) ) ;
      gsl_matrix_set( b , i , j , gsl_matrix_get( U , i + NROWS , j ) ) ;
    }
  }

 FREE :
  gsl_matrix_free( U ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( WORK ) ;

  return 0 ;
}

// set the a and b matrices
static void
set_ab( gsl_matrix *a , gsl_matrix *b ,
	const double *C0 , const double *C1 ,
	const size_t N , const size_t M )
{
  size_t i , j ;
  if( N == M ) {
    for( i = 0 ; i < N ; i++ ) {
      for( j = 0 ; j < N ; j++ ) {
	gsl_matrix_set( a , i , j , *C0 ) ; C0++ ;
	gsl_matrix_set( b , i , j , *C1 ) ; C1++ ;
      }
    }
  } else {
    svd_TLS( a , b , C0 , C1 , M , N ) ;
    /*
    size_t k ;
    // square it
    for( i = 0 ; i < N ; i++ ) {
      for( j = 0 ; j < N ; j++ ) {
	gsl_matrix_set( a , i , j , 0.0 ) ;
	gsl_matrix_set( b , i , j , 0.0 ) ;
      }
    }
    for( i = 0 ; i < N ; i++ ) {
      for( j = 0 ; j < N ; j++ ) {
	double l1 = 0.0 , l2 = 0.0 ;
	for( k = 0 ; k < M ; k++ ) {
	  l1 += C0[ i + k*N ] * C0[ j + k*N ] ;
	  l2 += C1[ i + k*N ] * C1[ j + k*N ] ;
	}
	gsl_matrix_set( a , i , j , l1 ) ;
	gsl_matrix_set( b , i , j , l2 ) ;

	//printf( " %f " , gsl_matrix_get( a , i , j ) ) ;
      }
      //printf( "\n" ) ;
    }
    */
    //
  }
  return ;
}

// solve the GEVP
struct resampled *
solve_GEVP( struct input_params *Input ,
	    const size_t N ,
	    const size_t M ,
	    const size_t t0 )
{
  if( N*M != Input -> Data.Nsim ) {
    fprintf( stderr , "[GEVP] N*M != %zu\n" , Input -> Data.Nsim ) ;
    return NULL ;
  }
  if( N > M ) {
    fprintf( stderr , "[GEVP] cannot solve when N states "
	     "are less than M correlators\n" ) ;
    return NULL ;
  }
  
  // initialise the generalised eigenvalues
  struct resampled *evalues = malloc( Input -> Data.Ndata[0] * N *
				      sizeof( struct resampled ) ) ;
  size_t i , j , k ;
  for( j = 0 ; j < Input ->Data.Ndata[0]*N ; j++ ) {
    evalues[j].resampled = malloc( Input -> Data.y[0].NSAMPLES *
				   sizeof( double ) ) ;
    evalues[j].restype = Input -> Data.y[0].restype ;
    evalues[j].NSAMPLES = Input -> Data.y[0].NSAMPLES ;
  }

  // set the data
  gsl_matrix *a  = gsl_matrix_alloc( N , N ) ;
  gsl_matrix *b  = gsl_matrix_alloc( N , N ) ;

  // allocations for GSL
  gsl_eigen_genv_workspace *work = gsl_eigen_genv_alloc( N ) ;
  gsl_vector_complex *alpha  = gsl_vector_complex_alloc( N ) ;
  gsl_vector *beta  = gsl_vector_alloc( N ) ;
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc( N , N ) ;
  
  // ugh loop order is all weird
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {

    double C0[ Input -> Data.Nsim ] , C1[ Input -> Data.Nsim ] ;

    fprintf( stdout , "[GEVP] solving GEVP %zu \n" , j ) ;
    
    for( k = 0 ; k < Input -> Data.y[0].NSAMPLES ; k++ ) {

      // put these in the linearised matrices
      size_t shift = 0 ;
      for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
	C0[i] = Input -> Data.y[ j  + shift ].resampled[k] ;
	C1[i] = Input -> Data.y[ t0 + shift ].resampled[k] ;
	shift += Input -> Data.Ndata[i] ;
      }

      // set AB by squaring the matrices, possibly
      set_ab( a , b , C0 , C1 , N , M ) ;
      
      // compute eigenvalues
      if( gevp( a , b , work , alpha , beta , evec ,
		N , j<t0 , false ) == FAILURE ) {
	fprintf( stderr , "[GEVP] GEVP solve failed \n" ) ;
        return NULL ;
      }

      // poke into solution
      for( i = 0 ; i < N ; i++ ) {
	evalues[j+Input->Data.Ndata[0]*i].resampled[ k ] =
	  ( gsl_vector_complex_get( alpha , i ).dat[0] )
	  / gsl_vector_get( beta , i ) ;
	//evalues[j+Input->Data.Ndata[0]*i].resampled[ k ] *= evalues[j+Input->Data.Ndata[0]*i].resampled[ k ] ;
      }
      //
    }

    // redo for the average
    size_t shift = 0 ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      C0[i] = Input -> Data.y[ j  + shift ].avg ;
      C1[i] = Input -> Data.y[ t0 + shift ].avg ;
      shift += Input -> Data.Ndata[i] ;
    }

    // set AB by squaring the matrices
    set_ab( a , b , C0 , C1 , N , M ) ;

    // compute eigenvalues
    if( gevp( a , b , work , alpha , beta , evec ,
	      N , j<t0 , true ) == FAILURE ) {
      fprintf( stderr , "[GEVP] GEVP solve failed \n" ) ;
      return NULL ;
    }

    // poke into solution
    for( i = 0 ; i < N ; i++ ) {
      evalues[j+Input->Data.Ndata[0]*i].avg =
	( gsl_vector_complex_get( alpha , i ).dat[0] )
	/ gsl_vector_get( beta , i ) ;
      //evalues[j+Input->Data.Ndata[0]*i].avg *= evalues[j+Input->Data.Ndata[0]*i].avg ;
    }
    //
  }

  gsl_matrix_complex_free( evec ) ;
  
  gsl_vector_complex_free( alpha ) ;
  gsl_vector_free( beta ) ;
  gsl_eigen_genv_free( work ) ;

  gsl_matrix_free( a ) ;
  gsl_matrix_free( b ) ;
  
  return evalues ;
}
