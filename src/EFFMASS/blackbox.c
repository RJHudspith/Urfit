/**
   @file blackbox.c
   @brief matrix-prony method for computing the effective mass

   From the paper what lattice theorists can learn from fMRI by Fleming
 */
#include <complex.h>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "gens.h"

// standard Prony's method
static void
blackbox1( double complex *x ,
	   const double *data ,
	   const size_t Ndata ,
	   const size_t Nstates ,
	   const size_t t ,
	   const size_t Nt )
{
  size_t i , j ;

  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( Nt , Nstates ) ;
  gsl_matrix *Q     = gsl_matrix_alloc( Nstates , Nstates ) ;
  gsl_vector *S     = gsl_vector_alloc( Nstates ) ;
  gsl_vector *Work  = gsl_vector_alloc( Nstates ) ;
  
  gsl_vector *beta  = gsl_vector_alloc( Nt ) ;
  gsl_vector *delta = gsl_vector_alloc( Nstates ) ;

  
  for( i = 0 ; i < Nt ; i++ ) {
    for( j = 0 ; j < Nstates ; j++ ) {
      const size_t didx = ( t + i + j )%Ndata ;
      gsl_matrix_set( alpha , i , j , -data[ didx ] ) ;
    }
  }

  for( i = 0 ; i < Nt ; i++ ) {
    const size_t yidx = ( t + Nstates + i )%Ndata ;
    gsl_vector_set( beta , i , data[ yidx ] ) ;
  }
  
  gsl_linalg_SV_decomp( alpha , Q , S , Work ) ;

  gsl_linalg_SV_solve( alpha , Q , S , beta , delta ) ;

  //printf( "States \n" ) ;
  double coeffs[ Nstates+1 ] ;

  for( i = 0 ; i < Nstates ; i++ ) {
    coeffs[ i ] = gsl_vector_get( delta , i ) ;
  }
  coeffs[Nstates] = 1 ;

  gsl_poly_complex_workspace *w =
    gsl_poly_complex_workspace_alloc( Nstates+1 ) ;
  double *z = malloc( 2 * Nstates * sizeof( double ) ) ;
  
  gsl_poly_complex_solve( coeffs , Nstates+1 , w , z ) ;
  
  //
  for( i = 0 ; i < Nstates ; i++ ) {
    x[ i ] = ( z[ 2*i ] + I * z[ 2*i+1 ] ) ;
  }

  gsl_poly_complex_workspace_free( w ) ;

  // free the gsl vectors
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( Work ) ;
  
  gsl_vector_free( beta ) ;
  gsl_vector_free( delta ) ;

  return ;
}

// standard Prony's method
static void
blackbox1p5( double complex *x ,
	     const double *data ,
	     const size_t Ndata ,
	     const size_t Nstates ,
	     const size_t t ,
	     const size_t Nt )
{
  size_t i , j ;

  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( Nt , Nt ) ;
  gsl_matrix *Q     = gsl_matrix_alloc( Nt , Nt ) ;
  gsl_vector *S     = gsl_vector_alloc( Nt ) ;
  gsl_vector *Work  = gsl_vector_alloc( Nt ) ;
  
  gsl_vector *beta  = gsl_vector_alloc( Nt ) ;
  gsl_vector *delta = gsl_vector_alloc( Nt ) ;

  
  for( i = 0 ; i < Nt ; i++ ) {
    for( j = 0 ; j < Nt ; j++ ) {
      const size_t didx = ( t + i + j )%Ndata ;
      gsl_matrix_set( alpha , i , j , -data[ didx ] ) ;
    }
  }

  for( i = 0 ; i < Nt ; i++ ) {
    const size_t yidx = ( t + Nt + i )%Ndata ;
    gsl_vector_set( beta , i , data[ yidx ] ) ;
  }
  
  gsl_linalg_SV_decomp( alpha , Q , S , Work ) ;

  for( i = Nstates ; i < Nt ; i++ ) {
    gsl_vector_set( S , i , 0.0 ) ;
  }
  
  gsl_linalg_SV_solve( alpha , Q , S , beta , delta ) ;

  //printf( "States \n" ) ;
  double coeffs[ Nt+1 ] ;

  for( i = 0 ; i < Nt ; i++ ) {
    coeffs[ i ] = gsl_vector_get( delta , i ) ;
  }
  coeffs[Nt] = 1 ;

  gsl_poly_complex_workspace *w =
    gsl_poly_complex_workspace_alloc( Nt+1 ) ;
  double *z = malloc( 2 * Nt * sizeof( double ) ) ;
  
  gsl_poly_complex_solve( coeffs , Nt+1 , w , z ) ;
  
  //
  for( i = 0 ; i < Nstates ; i++ ) {
    x[ i ] = ( z[ 2*(Nt-Nstates+i) ] + I * z[ 2*(Nt-Nstates+i)+1 ] ) ;
    //printf( "%zu %e \n" , i , creal( -clog( z[ 2*i ] ) ) ) ;
  }

  gsl_poly_complex_workspace_free( w ) ;

  // free the gsl vectors
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( Work ) ;
    
  gsl_vector_free( beta ) ;
  gsl_vector_free( delta ) ;

  return ;
}

// TLS variant
static void
blackbox3( double complex *x ,
	   const double *data ,
	   const size_t Ndata ,
	   const size_t Nstates ,
	   const size_t t ,
	   size_t Nt )
{
  size_t i , j ;

  if( Nt <= Nstates ) Nt += 1 ;

  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( Nt , Nstates+1 ) ;
  gsl_matrix *Q     = gsl_matrix_alloc( Nstates+1 , Nstates+1 ) ;
  gsl_vector *S     = gsl_vector_alloc( Nstates+1 ) ;
  gsl_vector *Work  = gsl_vector_alloc( Nstates+1 ) ;

  for( i = 0 ; i < Nt ; i++ ) {
    for( j = 0 ; j < Nstates ; j++ ) {
      const size_t didx = ( t + i + j )%Ndata ;
      gsl_matrix_set( alpha , i , j , -data[ didx ] ) ;
    }
    const size_t yidx = ( t + Nstates + i )%Ndata ;
    gsl_matrix_set( alpha , i , j , data[ yidx ] ) ;
  }
    
  gsl_linalg_SV_decomp( alpha , Q , S , Work ) ;

  const double NORM = -1.0 / gsl_matrix_get( Q , Nstates , Nstates ) ;
  if( fabs( NORM ) < 1E-15 ) {
    printf( "Singular NORM %e \n" , NORM ) ;
  }
  
  double coeffs[ Nstates+1 ] ;
  for( i = 0 ; i < Nstates ; i++ ) {
    coeffs[ i ] = NORM * gsl_matrix_get( Q , i , Nstates ) ;
  }
  coeffs[Nstates] = 1 ;
  
  // free the gsl vectors and matrices
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( Work ) ;
  
  gsl_poly_complex_workspace *w =
    gsl_poly_complex_workspace_alloc( Nstates+1 ) ;
  double *z = malloc( 2 * Nstates * sizeof( double ) ) ;
  
  if( gsl_poly_complex_solve( coeffs , Nstates+1 , w , z ) != GSL_SUCCESS ) {
    return ;
  }
  
  //
  for( i = 0 ; i < Nstates ; i++ ) {
    x[ i ] = z[ 2*i ] + I * z[ 2*i+1 ] ;
  }

  gsl_poly_complex_workspace_free( w ) ;

  return ;
}

// HTLS variant
static void
blackbox4( double complex *x ,
	   const double *data ,
	   const size_t Ndata ,
	   const size_t Nstates ,
	   const size_t t ,
	   const size_t Nt )
{
  size_t i , j ;

  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( Nt , Nstates+1 ) ;
  gsl_matrix *Q     = gsl_matrix_alloc( Nstates+1 , Nstates+1 ) ;
  gsl_vector *S     = gsl_vector_alloc( Nstates+1 ) ;
  gsl_vector *Work  = gsl_vector_alloc( Nstates+1 ) ;

  for( i = 0 ; i < Nt ; i++ ) {
    for( j = 0 ; j < Nstates ; j++ ) {
      const size_t didx = ( t + i + j )%Ndata ;
      gsl_matrix_set( alpha , i , j , -data[ didx ] ) ;
    }
    const size_t yidx = ( t + Nstates + i )%Ndata ;
    gsl_matrix_set( alpha , i , j , data[ yidx ] ) ;
  }

  gsl_linalg_SV_decomp( alpha , Q , S , Work ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( Work ) ;

  // truncate Hankel matrix to rank K, i.e. Nstates
  gsl_matrix *HK = gsl_matrix_alloc( 2*Nstates , Nstates ) ;
  for( i = 0 ; i < Nstates ; i++ ) {
    for( j = 0 ; j < Nstates ; j++ ) {
      gsl_matrix_set( HK , i , j ,
		      gsl_matrix_get( alpha , i , j ) ) ;
      gsl_matrix_set( HK , Nstates + i , j ,
		      gsl_matrix_get( alpha , i + 1 , j ) ) ;
    }
  }

  // free the gsl vectors
  gsl_matrix_free( alpha ) ;

  gsl_matrix *Qq     = gsl_matrix_alloc( Nstates , Nstates ) ;
  gsl_vector *Ss     = gsl_vector_alloc( Nstates ) ;
  gsl_vector *WorkQ  = gsl_vector_alloc( Nstates ) ;

  gsl_linalg_SV_decomp( HK , Qq , Ss , WorkQ ) ;

  gsl_matrix_free( Qq ) ;
  gsl_vector_free( Ss ) ;
  gsl_vector_free( WorkQ ) ;

  // split HK into two
  gsl_matrix *A = gsl_matrix_alloc( Nstates , Nstates ) ;
  gsl_matrix *B = gsl_matrix_alloc( Nstates , Nstates ) ;
  gsl_matrix *Binv = gsl_matrix_alloc( Nstates , Nstates ) ;
  gsl_matrix *C = gsl_matrix_alloc( Nstates , Nstates ) ;
  
  for( i = 0 ; i < Nstates ; i++ ) {
    for( j = 0 ; j < Nstates ; j++ ) {
      gsl_matrix_set( A , i , j , gsl_matrix_get( HK , i , j ) ) ;
      gsl_matrix_set( B , i , j , gsl_matrix_get( HK , i + Nstates , j  ) ) ;
    }
  }

  // invert B and multiply by A
  int signum ;
  gsl_permutation *perm = gsl_permutation_alloc( Nstates ) ;
  
  gsl_linalg_LU_decomp( B , perm , &signum ) ;
  gsl_linalg_LU_invert( B , perm , Binv ) ;

  for( i = 0 ; i < Nstates ; i++ ) {
    for( j = 0 ; j < Nstates ; j++ ) {
      register double loc_sum = 0.0 ;
      size_t k ;
      for( k = 0 ; k < Nstates ; k++ ) {
	loc_sum +=
	  gsl_matrix_get( A , i , k ) * gsl_matrix_get( Binv , k , j ) ;
      }
      gsl_matrix_set( C , i , j , loc_sum ) ;
    }
  }

  // compute the eigenvalues of C
  gsl_vector_complex *eval = gsl_vector_complex_alloc( Nstates ) ;
  gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc( Nstates ) ;
  gsl_eigen_nonsymm( C , eval , w ) ;

  for( i = 0 ; i < Nstates ; i++ ) {
    x[i] = GSL_REAL( gsl_vector_complex_get( eval , i ) ) + I * GSL_IMAG( gsl_vector_complex_get( eval , i ) ) ;
  }

  gsl_matrix_free( A ) ;
  gsl_matrix_free( B ) ;
  gsl_matrix_free( Binv  ) ;
  gsl_matrix_free( C  ) ;
  gsl_vector_complex_free( eval ) ;
  gsl_eigen_nonsymm_free( w ) ;

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

// qsort comparison
static int 
comp_sort( const void *elem1 , 
	   const void *elem2 ) 
{
  const double complex f = *( (double complex*)elem1 ) ;
  const double complex s = *( (double complex*)elem2 ) ;
  if( fabs(cimag(f)) > fabs(cimag(s)) ) return  1 ;
  if( fabs(cimag(f)) < fabs(cimag(s)) ) return -1 ;
  return 0 ;
}

/*
  So the data is very badly behaved, always
  what we want to do is have some data we trust and 
  some that we are unsure of

  worst offenders ( Nbad ) : solutions which are negative
                             or will yield a negative mass

  second worst ( Niffy ) : ones with imaginary parts
  
  we sort the real solutions into the trust list

  then we sort the imaginary ones

  and finally the worst ones
 */
static void
get_safe_masses( const size_t NSTATES ,
		 const size_t NDATA ,
		 double masses[ NSTATES ][ NDATA ] ,
		 const double complex *x ,
		 const size_t t )
{
  size_t Ntrust = 0 , Niffy = 0 , Nbad = 0 , i ;
  for( i = 0 ; i < NSTATES ; i++ ) {
    printf( "PRON_%zu %zu  %e %e \n" , i , t , creal( x[i] ) , cimag( x[i] ) ) ;
    if( creal( x[i] ) < 1.03 ) {
      Nbad++ ;
    } else if( fabs( cimag( x[i] ) ) < 1E-12 ) {
      Ntrust++ ;
    }
  }
  Niffy = NSTATES - Ntrust - Nbad ;

  printf( "Ntrust %zu , Niffy %zu , Nbad %zu \n" , Ntrust , Niffy , Nbad ) ;

  size_t t_idx = 0 , i_idx = 0 , b_idx = 0 ;
  double etrust[ Ntrust ] , eiffy[ Niffy ] , ebad[ Nbad ] ;
  for( i = 0 ; i < NSTATES ; i++ ) {
    if( creal( x[i] ) < 1.03 ) {
      ebad[ b_idx ] = fabs( creal( clog( x[i] ) ) ) ; b_idx++ ;
    } else if( fabs( cimag( x[i] ) ) < 1E-12 ) {
      etrust[ t_idx ] = fabs( creal( clog( x[i] ) ) ) ; t_idx++ ;
    } else {
      eiffy[ i_idx ] = fabs( creal( clog( x[i] ) ) ) ; i_idx++ ;
    }
  }
  printf( "Sols did \n" ) ;

  // sort them
  qsort( etrust , Ntrust , sizeof(double) , comp ) ;
  qsort( eiffy , Niffy, sizeof(double) , comp ) ;
  qsort( ebad , Nbad , sizeof(double) , comp ) ;

  size_t idx = 0 ;
  for( i = 0 ; i < Ntrust ; i++ ) {
    masses[idx][t] = etrust[i] ;
    idx++ ;
  }
  for( i = 0 ; i < Niffy ; i++ ) {
    masses[idx][t] = eiffy[i] ;
    idx++ ;
  }
  for( i = 0 ; i < Nbad ; i++ ) {
    masses[idx][t] = ebad[i] ;
    idx++ ;
  }

  for( i = 0 ; i < NSTATES ; i++ ) {
    printf( "Masses_%zu %e \n" , i , masses[i][t] ) ;
  }
  
  return ;
}

// compute the blackbox effective mass for a correlator of y-data
void
blackbox( const double *data ,
	  const size_t NDATA ,
	  const size_t NSTATES ,
	  double masses[ NSTATES ][ NDATA ] )
{ 
  // some cut off t
  const size_t smallt = NDATA ;

  double complex *x = malloc( NSTATES * sizeof( double complex ) ) ;

  size_t t ;
  for( t = 0 ; t < smallt ; t++ ) {

    blackbox1( x , data , NDATA , NSTATES , t , NSTATES+10 ) ;
    
    get_safe_masses( NSTATES , NDATA , masses , x , t ) ;
  }
  free( x ) ;

  return ;
}
