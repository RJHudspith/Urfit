/**
   @file svd.c
   @brief wrapper for gsl's SVD
 */
#include "gens.h"

//#define VERBOSE
#define SVD_COL_BALANCE
#define FAILURE 0
#define SUCCESS !FAILURE

#ifdef VERBOSE
// little print utility
static void
printmatrix( const double **mat ,
	     const size_t NROWS , 
	     const size_t NCOLS )
{
  printf( "\n" ) ;
  int i , j ;
  for( i = 0 ; i < NROWS ; i++ ) {
    printf( "|" ) ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      printf( " %f " , mat[i][j] ) ;
    }
    printf( "|\n" ) ;
  }
  printf( "\n" ) ;
  return ;
}
#endif

int 
svd_inverse( double **Ainv , 
	     const double **A ,
	     const size_t NCOLS ,
	     const size_t NROWS )
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
  register double sum , diff , tmp ;
  size_t i , j , k , FLAG = SUCCESS ;

  for( i = 0 ; i < NROWS ; i++ ){
    for( j = 0 ; j < NCOLS ; j++ ) {
      gsl_matrix_set( m , i , j , A[i][j] ) ;
    }
  }

  // do the decomposition
#ifdef SVD_COL_BALANCE
  gsl_vector *D = gsl_vector_alloc( NCOLS ) ;
  gsl_linalg_balance_columns( m , D ) ;
#endif

  if( gsl_linalg_SV_decomp( m , Q , S , WORK ) ) {
    printf( "GSL SVD comp failure \n" ) ;
    FLAG = FAILURE ;
    goto FREE ;
  }

  // loop these set 1.0 / Diag[i] to 0.0 if Diag[i] is crazy small
  for( i = 0 ; i < NCOLS ; i++ ) {
    Diag[ i ] = gsl_vector_get( S , i ) ; 
    tmp = 1.0 / Diag[i] ;
    InvDiag[ i ] = fabs( Diag[i] ) < 1.E-32 ? 0.0 : tmp ;
    #ifdef VERBOSE
    printf( "SVD %d %le\n" , i , Diag[i] ) ;
    printf( "SVD %d %le\n" , i , InvDiag[i] ) ;
    #endif
  }

  // test the solution to make sure it isn't too bad
  for( i = 0 ; i < NROWS ; i++ ) {
    for( k = 0 ; k < NCOLS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	sum += gsl_matrix_get( m , i , j ) * Diag[j] * gsl_matrix_get( Q , k , j ) ;
      }
      #ifdef SVD_COL_BALANCE
      diff += fabs( sum - A[i][k] / gsl_vector_get( D , k ) ) ;
      #else
      diff += fabs( sum - A[i][k] ) ;
      #endif
    }
  }

  // tell us how good the solution is and cry if it is too bad
#ifdef VERBOSE
  printf( "Decomposition accuracy :: %le \n" , diff ) ;
#endif
  diff /= (double)( NCOLS * NROWS ) ;
  if( diff > 1E-8 ) {
    printf( "SVD accuracy considered too low %e \n" , diff ) ;
    FLAG = FAILURE ;
    goto FREE ;
  } 

  // compute the product Ainv = V * ( 1.0 / Diag ) * U^T
  for( i = 0 ; i < NCOLS ; i++ ) {
    for( k = 0 ; k < NROWS ; k++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < NCOLS ;j++ ) {
	if( InvDiag[ j ] != 0.0 ) {
	  sum += gsl_matrix_get( Q , i , j ) * InvDiag[j] * gsl_matrix_get( m , k , j ) ;
	}
      }
      #ifdef SVD_COL_BALANCE
      Ainv[i][k] = sum / gsl_vector_get( D , i ) ;
      #else
      Ainv[i][k] = sum ;
      #endif
    }
  }

#ifdef VERBOSE
  // test our solution, should equal the NxN identity
  double **sol = malloc( NCOLS * sizeof( double* ) ) ;
  for( i = 0 ; i < NCOLS ; i++ ) {
    sol[i] = (double*)malloc( NCOLS * sizeof( double ) ) ;
    for( j = 0 ; j < NCOLS ; j++ ) {
      sum = 0.0 ;
      for( k = 0 ; k < NROWS ; k++ ) {
	sum += Ainv[i][k] * A[k][j] ;
      }
      sol[i][j] = sum ;
    }
  }
  printmatrix( (const double**)A , NROWS , NCOLS ) ;
  printmatrix( (const double**)Ainv , NCOLS , NROWS ) ;
  printmatrix( (const double**)sol , NCOLS , NCOLS ) ;
  for( i = 0 ; i < NCOLS ; i++ ) {
    free( sol[i] ) ;
  }
  free( sol ) ;
#endif

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
