/**
   @file correlation.c
   @brief compute the correlation matrix
 */
#include "gens.h"

#include "correlation.h"
#include "fitfunc.h"
#include "svd.h"

// fill the lower triangular part
static void
fill_lower_triangular( double **correlation ,
		       const size_t NDATA )
{
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < NDATA ; i++ ) {
    size_t j ;
    for( j = 0 ; j < i ; j++ ) {
      correlation[i][j] = correlation[j][i] ;
    }
  }
  return ;
}

// computes the upper section
static void
compute_upper_correlation( double **correlation , 
			   const struct resampled *data ,
			   const size_t NDATA ,
			   const corrtype CORRFIT )
{
  const size_t NSAMPLES = data[0].NSAMPLES ;

  // compute the correct normalisation
  double NORM = 1.0 ;
  size_t i ;

  // get the norm correct
  switch( data[0].restype ) {
  case RAWDATA :
    NORM = 1.0 / (double)( NSAMPLES * ( NSAMPLES - 1.0 ) ) ;
    break ;
  case JACKDATA :
    NORM = ( NSAMPLES - 1.0 ) / (double)NSAMPLES ;
    break ;
  case BOOTDATA :
    NORM = 1.0 / (double)NSAMPLES ;
    break ;
  }

  // set to zero
#pragma omp parallel for private(i)
  for( i = 0 ; i < NDATA ; i++ ) {
    size_t j ;
    for( j = 0 ; j < NDATA ; j++ ) {
      correlation[i][j] = 0.0 ;
    }
  }

  // diagonal one is pretty straightforward
  switch( CORRFIT ) {
  case UNWEIGHTED :
    break ;
  case UNCORRELATED :
#pragma omp parallel for private(i)
    for( i = 0 ; i < NDATA ; i++ ) {     
      const register double ave = data[i].avg ;
      register double sum = 0.0 ;
      size_t k ;
      for( k = 0 ; k < NSAMPLES ; k++ ) {
	sum += 
	  ( data[i].resampled[k] - ave ) *  
	  ( data[i].resampled[k] - ave ) ;
      }
      correlation[i][i] = sum * NORM ;
    }
    break ;
  case CORRELATED :
#pragma omp parallel for private(i)
    for( i = 0 ; i < NDATA ; i++ ) {
      const double avei = data[i].avg ;
      size_t j , k ;
      for( j = i ; j < NDATA ; j++ ) {
	const double avej = data[j].avg ;
	register double sum = 0.0 ;
	for( k = 0 ; k < NSAMPLES ; k++ ) {
	  sum += 
	    ( data[i].resampled[k] - avei ) *  
	    ( data[j].resampled[k] - avej ) ;
	}
	correlation[i][j] = sum * NORM ;
      }
    }
    break ;
  }
  return ;
}

// compute the covariance matrix, normalised by sigma_i sigma_j
void
correlations( double **correlation , 
	      const struct resampled *data ,
	      const corrtype CORRFIT ,
	      const size_t NDATA )
{
  // compute this correlation matrix
  compute_upper_correlation( correlation , data , NDATA , CORRFIT ) ;

  // divide by sigma^2
  modified_covariance( correlation , NDATA ) ;

  // correlation matrix is symmetric
  fill_lower_triangular( correlation , NDATA ) ;

  // write the correlation matrix?
  write_corrmatrix( correlation , NDATA ) ;

  return ;
}

// compute the inverse of the correlation matrix
double **
correlations_inv( const struct resampled *data ,
		  const corrtype CORRFIT ,
		  const size_t NDATA )
{
  double **Cinv = malloc( NDATA * sizeof( double* ) ) ;
  double **C    = malloc( NDATA * sizeof( double* ) ) ;

  // set Cinv to the identity just in case
  size_t i , j ;
  for( i = 0 ; i < NDATA ; i++ ) {
    Cinv[ i ] = malloc( NDATA * sizeof( double ) ) ;
    C[ i ]    = malloc( NDATA * sizeof( double ) ) ;
    for( j = 0 ; j < NDATA ; j++ ) {
      if( j == i ) {
	Cinv[ i ][ j ] = C[ i ][ j ] = 1.0 ;
      } else {
	Cinv[ i ][ j ] = C[ i ][ j ] = 0.0 ;
      }
    }
  }

  // compute this correlation matrix (false means not diagonal)
  compute_upper_correlation( C , data , NDATA , CORRFIT ) ;

  // correlation matrix is symmetric
  fill_lower_triangular( C , NDATA ) ;

  switch( CORRFIT ) {
  case UNWEIGHTED : break ;
  case  UNCORRELATED :
    for( i = 0 ; i < NDATA ; i++ ) {
      Cinv[i][i] = 1.0 / C[i][i] ;
    }
    break ;
  case CORRELATED : // compute the full inverse
    svd_inverse( Cinv , (const double**)C , NDATA , NDATA ) ;
    break ;
  }

  // free the correlation matrix
  for( i = 0 ; i < NDATA ; i++ ) {
    free( C[i] ) ;
  }
  free( C ) ;

  return Cinv ;
}

void
write_corrmatrix( double **correlation ,
		  const size_t NCUT )
{
  size_t i , j ;
  for( i = 0 ; i < NCUT ; i++ ) {
    for( j = 0 ; j < NCUT ; j++ ) {
      printf( "%f " , correlation[i][j] ) ;
    }
    printf( "\n" ) ;
  }
  return ;
}

void
write_corrmatrix_mathematica( double **correlation ,
			      const size_t NCUT )
{
  size_t i , j ;
  printf( "\n{{" ) ;
  for( i = 0 ; i < NCUT-1 ; i++ ) {
    for( j = 0 ; j < NCUT-1 ; j++ ) {
      printf( "%1.10f," , correlation[i][j] ) ;
    }     
    printf( "%1.10f},{\n" , correlation[i][j] ) ;
  }
  for( j = 0 ; j < NCUT-1 ; j++ ) {
    printf( "%1.10f," , correlation[i][j] ) ;
  }
  printf( "%1.10f}}\n" , correlation[i][j] ) ;
  return ;
}

void
write_corrmatrix_to_file( FILE *outfile , 
			  const double **correlation ,
			  const size_t NCUT )
{
  size_t i , j ;
  for( i = 0 ; i < NCUT ; i++ ) {
    for( j = 0 ; j < NCUT ; j++ ) {
      fprintf( outfile , "%1.3E " , correlation[i][j] ) ;
    }
    fprintf( outfile , "\n" ) ;
  }
  return ;
}
