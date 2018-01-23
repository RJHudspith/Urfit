/**
   @file GLS.c
   @brief generalised least squares with Tikhonov regularised priors
 */
#include "gens.h"

#include "chisq.h"
#include "svd.h"

#define WITH_SVD

// perform a generalised least squares iteration
int
gls_iter( void *fdesc ,
	  const void *data ,
	  const double **W ,
	  const double TOL )
{
  // point to the fit function
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;

  // check the matrix is there
  if( Fit -> f.U == NULL ) {
    fprintf( stderr , "[GLS] U matrix not initialised \n" ) ;
    return FAILURE ;
  }
  
  // point at the data
  struct data *Data = (struct data*)data ;

  // set the priors
  Fit -> f.Prior = Fit -> Prior ;

  // number of poly coefficients is Nlogic
  const size_t M = Fit -> f.N ;
  const size_t N = Fit -> Nlogic ;

  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( N , N ) ;
  gsl_vector *beta  = gsl_vector_alloc( N ) ;
  gsl_vector *delta = gsl_vector_alloc( N ) ;
  gsl_permutation *perm = gsl_permutation_alloc( N ) ;

  // usual counters
  size_t i , j , k , l ;

  // get the matrix description of our x data
  Fit -> linmat( Fit -> f.U , data , Fit -> Nparam  , Fit -> Nlogic ) ;

#ifdef VERBOSE
  // tell us what it looks like
  for( i = 0 ; i < M ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      printf( " %f " , Fit -> f.U[i][j] ) ;
    }
    printf( "\n" ) ;
  }
#endif
  
  // initialise alpha as (U^T W^-1 U + Q ) where Q is the diagonal matrix of prior errors
  for( i = 0 ; i < N ; i++ ) {
    for( j = i ; j < N ; j++ ) {
      register double sum = 0.0 ;
      for( k = 0 ; k < M ; k++ ) {
	switch( Fit -> f.CORRFIT ) {
	case UNWEIGHTED :
	  sum += Fit -> f.U[k][i] * Fit -> f.U[k][j] ;
	  break ;
	case UNCORRELATED :
	  sum += Fit -> f.U[k][i] * W[0][k] * Fit -> f.U[k][j] ;
	  break ;
	case CORRELATED :
	  for( l = 0 ; l < M ; l++ ) {
	    sum += Fit -> f.U[k][i] * W[k][l] * Fit -> f.U[l][j] ;
	  }
	  break ;
	}
      }
      gsl_matrix_set( alpha , i , j , sum ) ;
	
      // add in the priors
      if( i == j ) {
	if( Fit -> f.Prior[i].Initialised == true ) {
	  const double alpha_ii = gsl_matrix_get( alpha , i , i ) ;
	  const double sigmasq  = Fit -> f.Prior[i].Err * Fit -> f.Prior[i].Err ;
	  gsl_matrix_set( alpha , i , i ,
			  ( 1.0 / sigmasq ) + alpha_ii ) ;
	}
      } else {
	gsl_matrix_set( alpha , j , i , sum ) ;
      }
      //
    }
  }

#ifdef VERBOSE
  printf( "\n" ) ;
  for( i = 0 ; i < N ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      printf( " %1.2f " , gsl_matrix_get( alpha , i , j ) ) ;
    }
    printf( "\n" ) ;
  }
#endif

  // compute beta == U^T W^-1 y with tikhonov regularisation for the priors
  for( i = 0 ; i < N ; i++ ) {
    register double sum = 0.0 ;
    for( j = 0 ; j < M ; j++ ) {
      switch( Fit -> f.CORRFIT  ) {
      case UNWEIGHTED :
	sum += Fit -> f.U[j][i] * Data -> y[j] ;
	break ;
      case UNCORRELATED :
	sum += Fit -> f.U[j][i] * W[0][j] * Data -> y[j] ;
	break ;
      case CORRELATED :
	for( l = 0 ; l < M ; l++ ) {
	  sum += Fit -> f.U[j][i] * W[j][l] * Data -> y[l] ;
	}
	break ;
      }
    }
    if( Fit -> f.Prior[i].Initialised == true ) {
      sum += Fit -> f.Prior[i].Val / ( Fit -> f.Prior[i].Err * Fit -> f.Prior[i].Err ) ;
    }
    gsl_vector_set( beta , i , sum ) ;
  }

#ifdef SVD
  double A[N][N] , Ainv[N][N] ;
  for( i = 0 ; i < N ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      A[i][j] = gsl_matrix_get( alpha , i , j ) ;
    }
  }
  
  svd_inverse( Ainv , A , N , N , 1E-16 , true ) ;

  for( i = 0 ; i < N ; i++ ) {
    register double sum = 0.0 ;
    for( j = 0 ; j < N ; j++ ) {
      sum += Ainv[i][j] * gsl_vector_get( beta , i ) ;
    }
    Fit -> f.fparams[i] = sum ;
  }
  
#else
  // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
  int signum , Flag = SUCCESS ;
  if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
    fprintf( stderr , "[POLY_COEFF] LU decomp failed, rolling back to SVD\n" ) ;
    Flag = FAILURE ;
  }
  // if that worked we do the LU solve
  if( Flag != FAILURE ) {
    if( gsl_linalg_LU_solve( alpha , perm , beta , delta ) != GSL_SUCCESS ) {
      fprintf( stderr , "[POLY_COEFF] LU solve failure \n" ) ;
      Flag = FAILURE ;
    } else {
      // set the coefficients
      for( i = 0 ; i < N ; i++ ) {
        Fit -> f.fparams[i] = gsl_vector_get( delta , i ) ;
      }
    }
  }
#endif

  // set the parameter "f" and compute the chisq
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

  // free the gsl vectors
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;
  
  return Flag ;
}
