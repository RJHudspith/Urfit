/**
   @file svd.c
   @brief wrapper for gsl's SVD
 */
#include "gens.h"

int 
svd_inverse( double **Ainv , 
	     const double **A ,
	     const size_t NCOLS ,
	     const size_t NROWS ,
	     const double tolerance ,
	     const bool col_balance )
{
  // initial case dies
  if( NROWS < NCOLS ) {
    printf( "[SVD] Inverse | Rows > Columns "
	    ":: %zu vs %zu \n" , NROWS , NCOLS ) ;
    return FAILURE ;
  }

  // allocations
  gsl_matrix *m  = gsl_matrix_alloc( NROWS , NCOLS ) ;
  gsl_matrix *Q  = gsl_matrix_alloc( NCOLS , NCOLS ) ;
  gsl_vector *S  = gsl_vector_alloc( NCOLS ) ;
  gsl_vector *WORK  = gsl_vector_alloc( NCOLS ) ;
  double Diag[ NCOLS ] , InvDiag[ NCOLS ] ;
  register double sum , diff = 0.0 , tmp ;
  size_t i , j , k , FLAG = SUCCESS ;

  // initialize our matrix
  for( i = 0 ; i < NROWS ; i++ ){
    for( j = 0 ; j < NCOLS ; j++ ) {
      gsl_matrix_set( m , i , j , A[i][j] ) ;
    }
  }

  // balance the columns if we want
  gsl_vector *D = NULL ;
  if( col_balance == true ) {
    D = gsl_vector_alloc( NCOLS ) ;
    gsl_linalg_balance_columns( m , D ) ;
  }

  // perform the decomposition
  if( gsl_linalg_SV_decomp( m , Q , S , WORK ) ) {
    printf( "[SVD] GSL SVD comp failure \n" ) ;
    FLAG = FAILURE ;
    goto FREE ;
  }

  // loop these set 1.0 / Diag[i] to 0.0 if Diag[i] is crazy small
  size_t nsmall = 0 ;
  for( i = 0 ; i < NCOLS ; i++ ) {
    Diag[ i ] = gsl_vector_get( S , i ) ;
    printf( "[SVD] diag %zu -> %e \n" , i , Diag[i] ) ; 
    if( Diag[i] < tolerance ) {
      tmp = 0.0 ;
      nsmall++ ;
    } else {
      tmp = 1.0 / Diag[i] ;
    }
    // make sure it doesn't get ridiculously big either
    InvDiag[ i ] = fabs( Diag[i] ) < 1E-40 ? 0.0 : tmp ;
  }
  fprintf( stdout , "[SVD] svd cut %1.2f percent omitted\n" ,
	   100*(nsmall)/(double)NCOLS ) ;
  if( nsmall == NCOLS ) {
    fprintf( stderr ,
	     "[SVD] cut too aggressive, covariance matrix singular\n" ) ;
    return FAILURE ;
  }

  // test the solution to make sure it isn't too bad
  diff = 0.0 ;
  for( i = 0 ; i < NROWS ; i++ ) {
    for( k = 0 ; k < NCOLS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	sum += gsl_matrix_get( m , i , j ) * Diag[j] * gsl_matrix_get( Q , k , j ) ;
      }
      // if we balanced the columns the shift is a little different
      if( col_balance == true ) {
	diff += fabs( sum - A[i][k] / gsl_vector_get( D , k ) ) ;
      } else {
	diff += fabs( sum - A[i][k] ) ;
      }
    }
  }

  // tell us how good the solution is and cry if it is too bad
#ifdef VERBOSE
  printf( "[SVD] Decomposition accuracy :: %le \n" , diff ) ;
#endif
  diff /= (double)( NCOLS * NROWS ) ;
  if( diff > 1E-14*Diag[0] ) {
    printf( "[SVD] accuracy below tolerance %e < %e \n" , diff*Diag[0] ,
	    1E-12 ) ;
    FLAG = FAILURE ;
    goto FREE ;
  } 

  // compute the product Ainv = V * ( 1.0 / Diag ) * U^T
  for( i = 0 ; i < NCOLS ; i++ ) {
    for( k = 0 ; k < NROWS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	sum += gsl_matrix_get( Q , i , j ) * InvDiag[j] * gsl_matrix_get( m , k , j ) ;
      }
      // rebalance the columns
      if( col_balance == true ) {
	Ainv[i][k] = sum / gsl_vector_get( D , i ) ;
      } else {
	Ainv[i][k] = sum ;
      }
    }
  }

 FREE :
  gsl_matrix_free( m ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( WORK ) ;

#ifdef SVD_COL_BALANCE
  gsl_vector_free( D ) ;
#endif

  return FLAG ;
}
