/**
   @file pade_laplace.c
   @brief laplace-pade method
 */
#include "gens.h"

#include <gsl/gsl_poly.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "pade_coefficients.h"

// compute the coefficients
static void
dLdp( double *d ,
      const double *fac ,
      const size_t Nders ,
      const double *x ,
      const double *y ,
      const size_t N ,
      const double p0 )
{
  double epx[ N ] ;
  size_t i , j ;
  // initialise the exponential factors and the powers of t
  for( i = 0 ; i < N ; i++ ) {
    epx[ i ] = exp( -p0 * x[i] ) * y[i] ;
  }
  for( i = 0 ; i < Nders ; i++ ) {
    d[i] = 0.0 ;
    for( j = 0 ; j < N-1 ; j++ ) {
      d[i] += ( epx[j+1] + epx[j] ) * ( x[j+1] - x[j] ) ;
      epx[j] *= -x[j] ;
    }
    epx[j] *= -x[j] ;
    d[i] /= ( 2 * fac[i] ) ;
    #ifdef VERBOSE
    printf( "D[ %zu ] %e %e \n" , i , d[i] , fac[i] ) ;
    #endif
  }
  return ;
}

// qsort comparison
static int 
comp( const void *elem1 , 
      const void *elem2 ) 
{
  const double f = *( (double*)elem1 ) ;
  const double s = *( (double*)elem2 ) ;
  if (f > s) return  1 ;
  if (f < s) return -1 ;
  return 0 ;
}

static int
stable_pade( double *sims ,
	     const double *d ,
	     const size_t Nexp ,
	     const double p0 )
{
  double results[ Nexp ][ Nexp ] ;
  size_t i , j ;
  for( j = 0 ; j < Nexp ; j++ ) {

    sims[j] = 0.0 ;
    
    for( i = 0 ; i < Nexp ; i++ ) {
      results[j][i] = 0.0 ;
    }
    
    const size_t N = Nexp*Nexp + j ;
    
    double *pade = malloc( 2*N * sizeof( double ) ) ;
    pades_from_poly( pade , d , N , N ) ;

    double *denom = malloc( ( N + 1 ) * sizeof( double ) ) ;
    denom[0] = 1 ;
    for( i = 0 ; i < N  ; i++ ) {
      denom[ i + 1 ] = pade[ N + i ] ;
    }

    gsl_poly_complex_workspace *w =		\
      gsl_poly_complex_workspace_alloc( N+1 ) ;
    double *z = malloc( 2 * ( N ) * sizeof( double ) ) ;
  
    gsl_poly_complex_solve( denom , N+1 , w , z ) ;

    // set the real data
    size_t idx = 0 ;
    for( i = 0 ; i < N ; i++ ) {
      //#ifdef VERBOSE
      printf( "RESULTS :: %zu %zu -> %f %f \n",
	      j , i , z[2*i] + p0 , z[2*i+1] + p0) ;
      //#endif
      if( fabs( z[2*i+1] ) < 1E-6 && // is purely real
	  fabs( z[2*i] + p0 ) < 3 && // is small
	  ( z[2*i] + p0 ) < 0        // is less than zero
	  ) {
	if( idx < Nexp ) {
	  results[j][idx] = z[2*i] + p0 ;
	  #ifdef VERBOSE
	  printf( "Results post filter : %f \n" , results[j][idx] ) ;
	  #endif
	}
	idx++ ;
      }
    }
    gsl_poly_complex_workspace_free( w ) ;
    free( pade ) ;
    free( denom ) ;
    free( z ) ;
  }

  // quicksort the data
  int Flag = SUCCESS ;
  for( i = 0 ; i < Nexp ; i++ ) {
    qsort( results[i] , Nexp , sizeof(double) , comp ) ;
  }
  for( i = 0 ; i < Nexp ; i++ ) {
    size_t Ntot = 0 ;
    sims[i] = 0.0 ;
    for( j = 0 ; j < Nexp ; j++ ) {
      if( results[j][i] != 0.0 ) {
	sims[i] += results[j][i] ;
	Ntot++ ;
      }
    }
    if( Ntot == 0 ) {
      fprintf( stderr , "[PADE-LAPLACE] could not determine all of the pole"
	       " masses Problem with -> Exp %zu \n" , i ) ;
      Flag = FAILURE ;
    } else {
      sims[i] /= Ntot ;
    }
  }
  return Flag ;
}

// get amplitudes
static int
get_amplitudes( double *Amps ,
		const double *x ,
		const double *y ,
		const size_t Ndata ,
		const double *masses ,
		const double Nexps )
{
  // compute the amplitudes using
  // ( A^T A ) A^T y = x
  // where the matrix A is e^mx
  gsl_matrix *alpha = gsl_matrix_alloc( Nexps , Nexps ) ;
  gsl_vector *beta  = gsl_vector_alloc( Nexps ) ;
  gsl_vector *delta = gsl_vector_alloc( Nexps ) ;
  gsl_permutation *perm = gsl_permutation_alloc( Nexps ) ;

  size_t i , j , k ;
  for( i = 0 ; i < Nexps ; i++ ) {
    // set alpha
    register double sum ;
    for( j = i ; j < Nexps ; j++ ) {
      sum = 0.0 ;
      for( k = 0 ; k < Ndata ; k++ ) {
	sum += exp( ( masses[i] + masses[j] ) * x[k] ) ;
      }
      gsl_matrix_set( alpha , i , j , sum ) ;
      gsl_matrix_set( alpha , j , i , sum ) ;
    }
    // set beta
    sum = 0.0 ;
    for( k = 0 ; k < Ndata ; k++ ) {
      sum += y[k] * exp( masses[i] * x[k] ) ;
    }
    gsl_vector_set( beta , i , sum ) ;
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

  // tell us the result
  for( i = 0 ; i < Nexps ; i++ ) {
    fprintf( stdout , "%f e ^ { %f * x } \n" ,
	     gsl_vector_get( delta , i ) , masses[i] ) ;
    Amps[i] = gsl_vector_get( delta , i ) ;
  }

  // free the gsl vectors
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  return Flag ;
}

int
pade_laplace( double *fparams ,
	      const double *x ,
	      const double *y ,
	      const size_t Ndata ,
	      const size_t Nexps )
{
  const size_t Nders = 30 ;
  
  double *d = malloc( Nders * sizeof( double ) ) ;
  double *fac = malloc( Nders * sizeof( double ) ) ;
  double Masses[ Nexps ] , Amps[ Nexps ] ;
  const double p0 = 1.0 ;
  
  size_t i ;
  fac[0] = 1 ;
  for( i = 1 ; i < Nders ; i++ ) {
    fac[i] = fac[i-1] * i ;
  }

  // compute the derivatives of the amplitudes
  dLdp( d , fac , Nders , x , y , Ndata , p0 ) ;
  
  // get the poles of the pade representation
  int Flag ;
  Flag = stable_pade( Masses , d , Nexps , p0 ) ;
  if( Flag == FAILURE ) goto free ;
  
  // compute the amplitudes from the masses
  Flag = get_amplitudes( Amps , x , y , Ndata , Masses , Nexps ) ;
  if( Flag == FAILURE ) goto free ;
  
  // set the fit parameters
  for( i = 0 ; i < Nexps ; i++ ) {
    fparams[ 2*i + 0 ] = Amps[ i ] ;
    fparams[ 2*i + 1 ] = Masses[ i ] ;
  }

 free :
  free( d ) ;
  free( fac ) ;

  return Flag ;
}
