/**
   @file correlation.c
   @brief compute the correlation matrix
 */
#include "gens.h"

#include "correlation.h"
#include "fitfunc.h"
#include "svd.h"

// computes the upper section
static void
compute_upper_correlation( double **correlation , 
			   const struct resampled *data ,
			   const size_t NDATA ,
			   const corrtype CORRFIT )
{
  // initialise number of samples
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
      correlation[0][i] = sum * NORM ;
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

// divide by sigma^2
static void
modified_covariance( double **correlation ,
		     const size_t Ndata )
{
  // compute sigma
  double sigma[ Ndata ] ;
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < Ndata ; i++ ) {
    sigma[ i ] = sqrt( correlation[i][i] ) ;
  }
  
  // rescale correlation matrix by variance
#pragma omp parallel for private(i)
  for( i = 0 ; i < Ndata ; i++ ) {
    size_t j ;
    for( j = i ; j < Ndata ; j++ ) {
      // edge case that I quite often hit
      if( sigma[i] == 0.0 || sigma[j] == 0.0 ) {
	if( j == i ) {
	  correlation[i][i] = 1.0 ;
	} else {
	  correlation[i][i] = 0.0 ;
	}
      // otherwise rescale
      } else {
	correlation[i][j] /= ( sigma[i] * sigma[j] ) ;
      }
      //
    }
  }
  return ;
}

// compute the covariance matrix, normalised by sigma_i sigma_j
void
correlations( double **correlation , 
	      const struct resampled *data ,
	      const corrtype CORRFIT ,
	      const size_t Ndata ,
	      const bool Divided )
{
  // compute this correlation matrix
  compute_upper_correlation( correlation , data , Ndata , CORRFIT ) ;

  // divide by sigma^2
  modified_covariance( correlation , Ndata ) ;

  // correlation matrix is symmetric
  if( Divided == true ) {
    fill_lower_triangular( correlation , Ndata ) ;
  }
  
  // write the correlation matrix?
  write_corrmatrix( (const double**)correlation , Ndata ) ;

  return ;
}

// compute the inverse of the correlation matrix
double **
correlations_inv( const struct resampled *data ,
		  const corrtype CORRFIT ,
		  const size_t Ndata ,
		  const double Tolerance ,
		  const bool Divided )
{
  double **Cinv = malloc( Ndata * sizeof( double* ) ) ;
  double **C    = malloc( Ndata * sizeof( double* ) ) ;

  // set Cinv to the identity just in case
  size_t i , j ;
  for( i = 0 ; i < Ndata ; i++ ) {
    Cinv[ i ] = malloc( Ndata * sizeof( double ) ) ;
    C[ i ]    = malloc( Ndata * sizeof( double ) ) ;
    for( j = 0 ; j < Ndata ; j++ ) {
      if( j == i ) {
	Cinv[ i ][ j ] = C[ i ][ j ] = 1.0 ;
      } else {
	Cinv[ i ][ j ] = C[ i ][ j ] = 0.0 ;
      }
    }
  }

  // compute this correlation matrix (false means not diagonal)
  compute_upper_correlation( C , data , Ndata , CORRFIT ) ;

  // divides elements by sigma_i sigma_j
  if( Divided == true ) {
    modified_covariance( C , Ndata ) ;
  }
  
  // correlation matrix is symmetric
  fill_lower_triangular( C , Ndata ) ;

  switch( CORRFIT ) {
  case UNWEIGHTED : break ;
  case  UNCORRELATED :
    for( i = 0 ; i < Ndata ; i++ ) {
      Cinv[i][i] = fabs( C[i][i] ) > 1.E-32 ? 1.0 / C[i][i] : 0.0 ;
    }
    break ;
  case CORRELATED : // compute the full inverse
    svd_inverse( Cinv , (const double**)C ,
		 Ndata , Ndata , Tolerance , false ) ;
    break ;
  }

  // free the correlation matrix
  for( i = 0 ; i < Ndata ; i++ ) {
    free( C[i] ) ;
  }
  free( C ) ;

  return Cinv ;
}

// compute the inverse of the correlation matrix
int
inverse_correlation( struct data_info *Data ,
		     const struct fit_info Fit )
{
  // space for the correlation matrix
  double **C = NULL ;
  size_t i ;
  int flag = SUCCESS ;

  // if we don't weight by the correlation matrix we don't allocate
  // any space for it
  switch( Fit.Corrfit ) {
  case UNWEIGHTED : return flag ;
  case UNCORRELATED :
    C = calloc( 1 , sizeof( double* ) ) ;
    C[ 0 ] = calloc( Data -> Ntot , sizeof( double ) ) ;
    Data -> Cov.W = calloc( 1 , sizeof( double* ) ) ;
    Data -> Cov.W[ 0 ] = calloc( Data -> Ntot , sizeof( double ) ) ;

    // compute this correlation matrix (false means not diagonal)
    compute_upper_correlation( C , Data -> y , Data -> Ntot , Fit.Corrfit ) ;

    for( i = 0 ; i < Data -> Ntot ; i++ ) {
      Data -> Cov.W[0][i] = fabs( C[0][i] ) > 1.E-32 ? 1.0 / C[0][i] : 1.E-32 ;
    }

    for( i = 0 ; i < 1 ; i++ ) {
      free( C[0] ) ;
    }
    free( C ) ;
    
    break ;
  case CORRELATED :
    // allocate the correlation matrices
    C = calloc( Data -> Ntot , sizeof( double* ) ) ;
    Data -> Cov.W = calloc( Data -> Ntot , sizeof( double* ) ) ;
    
    // set W to the identity just in case
    for( i = 0 ; i < Data -> Ntot ; i++ ) {
      Data -> Cov.W[ i ] = calloc( Data -> Ntot , sizeof( double ) ) ;
      C[ i ] = calloc( Data -> Ntot , sizeof( double ) ) ;
      Data -> Cov.W[ i ][ i ] = C[ i ][ i ] = 1.0 ;
    }

    // compute this correlation matrix (false means not diagonal)
    compute_upper_correlation( C , Data -> y , Data -> Ntot , Fit.Corrfit ) ;

    // divides elements by sigma_i sigma_j
    if( Data -> Cov.Divided_Covariance == true ) {
      modified_covariance( C , Data -> Ntot ) ;
    }
  
    // correlation matrix is symmetric
    fill_lower_triangular( C , Data -> Ntot ) ;

    // compute the inverse by svd
    if( svd_inverse( Data -> Cov.W , (const double**)C ,
		     Data -> Ntot , Data -> Ntot ,
		     Data -> Cov.Eigenvalue_Tol ,
		     Data -> Cov.Column_Balanced ) == FAILURE ) {
      flag = FAILURE ;
    }

    // free the correlation matrix
    for( i = 0 ; i < Data -> Ntot ; i++ ) {
      free( C[i] ) ;
    }
    free( C ) ;
    
    break ;
  }

  return flag ;
}

void
write_corrmatrix( const double **correlation ,
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
write_corrmatrix_mathematica( const double **correlation ,
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
