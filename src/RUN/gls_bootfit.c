/**
   @file GLS.c
   @brief generalised least squares
 */
#include "gens.h"

#include "chisq.h"
#include "resampled_ops.h"
#include "stats.h"
#include "svd.h"

// compute a particular coefficient
int
single_gls( double *coeffs ,
	    double *chisq ,
	    const struct data_info Data ,
	    const struct fit_info Fit ,
	    const size_t sample_idx ,
	    const bool is_average )
{
  // number of poly coefficients is Nlogic
  const size_t M = Data.Ntot ;
  const size_t N = Fit.Nlogic ;

  size_t i , j , k , l ;
  
  // matrix allocations
  double **U = malloc( M * sizeof( double* ) ) ;
  for( i = 0 ; i < M ; i++ ) {
    U[i] = malloc( N * sizeof( double ) ) ;
    for( j = 0 ; j < N ; j++ ) {
      U[i][j] = 0.0 ;
    }
  }

  gsl_matrix *alpha = gsl_matrix_alloc( N , N ) ;
  gsl_vector *beta  = gsl_vector_alloc( N ) ;
  gsl_vector *delta = gsl_vector_alloc( N ) ;
  gsl_permutation *perm = gsl_permutation_alloc( N ) ;
  
  // allocate the matrix with increasing powers of x
  register double sum , x0 ;
  for( i = 0 ; i < M ; i++ ) {
    x0 = 1 ;
    for( j = 0 ; j < Fit.Nparam ; j++ ) {
      U[i][ Fit.map[i].p[j] ] = x0 ;
      if( is_average == true ) {
	x0 *= Data.x[i].avg ;
      } else {
	x0 *= Data.x[i].resampled[sample_idx] ;
      }
    }
  }

  #ifdef VERBOSE
  // tell us what it looks like
  for( i = 0 ; i < M ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      printf( " %f " , U[i][j] ) ;
    }
    printf( "\n" ) ;
  }
  #endif
  
  // initialise alpha as (U^T W^-1 U)
  for( i = 0 ; i < N ; i++ ) {
    for( j = i ; j < N ; j++ ) {
      register double sum = 0.0 ;
      for( k = 0 ; k < M ; k++ ) {
	switch( Fit.Corrfit ) {
	case UNWEIGHTED :
	  sum += U[k][i] * U[k][j] ;
	  break ;
	case UNCORRELATED :
	  sum += U[k][i] * Data.Cov.W[0][k] * U[k][j] ;
	  break ;
	case CORRELATED :
	  for( l = 0 ; l < M ; l++ ) {
	    sum += U[k][i] * Data.Cov.W[k][l] * U[l][j] ;
	  }
	  break ;
	}
      }
      gsl_matrix_set( alpha , i , j , sum ) ;
      gsl_matrix_set( alpha , j , i , sum ) ;
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

  // compute beta == U^T W^-1 y
  for( i = 0 ; i < N ; i++ ) {
    register double sum = 0.0 ;
    for( j = 0 ; j < M ; j++ ) {
      switch( Fit.Corrfit ) {
      case UNWEIGHTED :
	if( is_average == true ) {
	  sum += U[j][i] * Data.y[j].avg ;
	} else {
	  sum += U[j][i] * Data.y[j].resampled[ sample_idx ] ;
	}
	break ;
      case UNCORRELATED :
	if( is_average == true ) {
	  sum += U[j][i] * Data.Cov.W[0][j] * Data.y[j].avg ;
	} else {
	  sum += U[j][i] * Data.Cov.W[0][j] *
	    Data.y[j].resampled[ sample_idx ] ;
	}
	break ;
      case CORRELATED :
	for( l = 0 ; l < M ; l++ ) {
	  if( is_average == true ) {
	    sum += U[j][i] * Data.Cov.W[j][l] *
	      Data.y[l].avg ;
	  } else {
	    sum += U[j][i] * Data.Cov.W[j][l] *
	      Data.y[l].resampled[sample_idx] ;
	  }
	}
	break ;
      }
    }
    gsl_vector_set( beta , i , sum ) ;
  }
  
  // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
  bool svd_flag = false ;
  int signum , Flag = SUCCESS ;
  if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
    fprintf( stderr , "[POLY_COEFF] LU decomp failed, rolling back to SVD\n" ) ;
    svd_flag = true ;
  }

  // if that worked we do the LU solve
  if( svd_flag != true ) {
    if( gsl_linalg_LU_solve( alpha , perm , beta , delta ) != GSL_SUCCESS ) {
      fprintf( stderr , "[POLY_COEFF] LU solve failure \n" ) ;
      svd_flag = true ;
    } else {
      // set the coefficients
      for( i = 0 ; i < N ; i++ ) {
	coeffs[i] = gsl_vector_get( delta , i ) ;
      }
    }
  }

  if( svd_flag == true ) {
    printf( "[POLY] regressing to the SVD \n" ) ;
    
    double **UtU_inv = malloc( N * sizeof( double* ) ) ;
    double **UtU    = malloc( N * sizeof( double* ) ) ;
    for( i = 0 ; i < N ; i++ ) {
      UtU_inv[i] = malloc( N * sizeof( double ) ) ;
      UtU[i]     = malloc( N * sizeof( double ) ) ;
      for( j = 0 ; j < N ; j++ ) {
	UtU[i][j] = gsl_matrix_get( alpha , i , j ) ;
      }
    }
    // if the SVD screws up we set failure and free allocations
    if( svd_inverse( UtU_inv , (const double**)UtU , N , N , 1E-14 , true
		     ) == FAILURE ) {
      Flag = FAILURE ;
    }
    
    // multiply the inverse by the data to obtain the coefficients
    for( i = 0 ; i < N ; i++ ) {
      sum = 0.0 ;
      for( j = 0 ; j < N ; j++ ) {
	sum += UtU_inv[i][j] * gsl_vector_get( beta , j ) ;
      }
      coeffs[i] = sum ;
    }

    // memfree
    for( i = 0 ; i < M ; i++ ) {
      free( UtU_inv[i] ) ;
      free( UtU[i] ) ;
    }
    free( UtU_inv ) ;
    free( UtU ) ;
  }

  *chisq = 0.0 ;
  // compute the full chisq here
  // chisq = f_i W_ij f_j
  // where f_i = p(x_i) - y_i
  // and p(x_i) is our poly solution at x

  double f[ M ] ;   // compute a little f
  for( i = 0 ; i < M ; i++ ) {
    register double sum = 0.0 ;
    for( j = 0 ; j < N ; j++ ) {
      sum += U[i][j] * coeffs[j] ;
    }
    if( is_average == true ) {
      f[i] = sum - Data.y[i].avg ;
    } else {
      f[i] = sum - Data.y[i].resampled[ sample_idx ] ;
    }
  }

  // compute the chisq from it
  *chisq = 0.0 ;
  for( i = 0 ; i < M ; i++ ) {
    switch( Fit.Corrfit ) {
    case UNWEIGHTED :
      *chisq += f[i] * f[i] ;
      break ;
    case UNCORRELATED :
      *chisq += f[i] * Data.Cov.W[0][i] * f[i] ;
      break ;
    case CORRELATED :
      for( j = 0 ; j < M ; j++ ) {
	*chisq += f[i] * Data.Cov.W[i][j] * f[j] ;
      }
      break ;
    }
  }

  // free the gsl vectors
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  // free this stuff
  for( i = 0 ; i < M ; i++ ) {
    free( U[i] ) ;
  }
  free( U ) ;
  
  return Flag ;
}

// computes poly coefficients using ( U W^-1 U^T ) ( U^T W^-1 y )
struct resampled *
GLS( const struct data_info Data ,
     const struct fit_info Fit )
{
  const size_t N = Fit.Nlogic ;

  // solve by LU? rolls back to column-balanced SVD if it can't
  size_t i , j ;

  struct resampled *fitparams = malloc( N * sizeof( struct resampled ) ) ;
  for( i = 0 ; i < N ; i++ ) {
    fitparams[i] = init_dist( NULL , Data.y[0].NSAMPLES ,
			      Data.y[0].restype ) ;
    fitparams[i].avg = UNINIT_FLAG ;
  }
  struct resampled chisq = init_dist( NULL , Data.y[0].NSAMPLES ,
				      Data.y[0].restype ) ;

  // do the average
  double ave_coeffs[ N ] , ave_chi = 0.0 ;
  single_gls( ave_coeffs , &ave_chi , Data , Fit , 0 , true ) ;
  for( j = 0 ; j < N ; j++ ) {
    fitparams[j].avg = ave_coeffs[j] ;
  }
  chisq.avg = ave_chi ;

  // do all of the individual bootstraps
#pragma omp parallel for private(i)
  for( i = 0 ; i < chisq.NSAMPLES ; i++ ) {
    double coeffs[ N ] , chi = 0.0 ;
    single_gls( coeffs , &chi , Data , Fit , i , false ) ;
    for( j = 0 ; j < N ; j++ ) {
      fitparams[j].resampled[i] = coeffs[j] ;
    }
    chisq.resampled[i] = chi ;
  }
  
  divide_constant( &chisq , ( Data.Ntot - Fit.Nlogic ) ) ;

  printf( "[CHISQ / (d.o.f)] %e %e \n" , chisq.avg , chisq.err ) ;

  // tell us what we have computed
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    compute_err( &fitparams[i] ) ;
    if( i < Fit.Nparam ) {
      if( Fit.Sims[i] == true ) {
	printf( "-> SIMUL " ) ;
      }
    }
    printf( "PARAM_%zu %f %f \n" , i , fitparams[i].avg , fitparams[i].err ) ;
  }

  free( chisq.resampled ) ;
  
  return fitparams ;
}
