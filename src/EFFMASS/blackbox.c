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

#include "gens.h"

// standard Prony's method
static void
blackbox2( double complex *x ,
	   const double *data ,
	   const size_t Ndata ,
	   const size_t Nstates ,
	   const size_t t ,
	   const size_t Nt )
{
  size_t i , j , k , n ;
  //int signum ;
  double *Y = malloc( Nt*sizeof( double ) ) ;
  double Y2[ Nt ][ Nstates ] ;
  for( i = 0 ; i < Nt ; i++ ) {
    const size_t yidx = ( t + Nstates + i )%Ndata ;
    
    for( j = 0 ; j < Nstates ; j++ ) {
      const size_t didx = ( t + i + j )%Ndata ;
      Y2[ i ][ j ] = -data[ didx ] ;
    }
    Y[i] = data[ yidx ] ;
  }

  // make it into a square Nstates*Nstates matrix by left multiplying by Y2^T
  // gsl matrix allocations
  gsl_matrix *alpha = gsl_matrix_alloc( Nstates , Nstates ) ;
  gsl_matrix *Q     = gsl_matrix_alloc( Nstates , Nstates ) ;
  gsl_vector *S     = gsl_vector_alloc( Nstates ) ;
  gsl_vector *Work  = gsl_vector_alloc( Nstates ) ;
  
  gsl_vector *beta  = gsl_vector_alloc( Nstates ) ;
  gsl_vector *delta = gsl_vector_alloc( Nstates ) ;
  
  for( i = 0 ; i < Nstates ; i++ ) {
    register double sum = 0.0 ;
    for( j = 0 ; j < Nstates ; j++ ) {
      sum = 0.0 ;
      for( k = 0 ; k < Nt ; k++ ) {
	sum += Y2[ k ][ i ] * Y2[ k ][ j ] ;
      }
      gsl_matrix_set( alpha , i , j , sum ) ;
    }
    sum = 0.0 ;
    for( k = 0 ; k < Nt ; k++ ) {
      sum += Y2[ k ][ i ] * Y[k] ;
    }
    gsl_vector_set( beta , i , sum ) ;
  }

#if 0
  printf( "%f \n" , data[t] ) ;
  printf( "%f \n" , data[(t+1)%Ndata] ) ;

  // have a look at it
  printf( "Ymatrix\n" ) ;
  for( i = 0 ; i < Nstates ; i++ ) {
    for( j = 0 ; j < Nstates ; j++ ) {
      printf( " %f " , gsl_matrix_get( alpha , i , j ) ) ;
    }
    printf( "\n" ) ;
  }
  
  printf( "Ys\n" ) ;
  for( i = 0 ; i < Nstates ; i++ ) {
    printf( " %f " , gsl_vector_get( beta , i ) ) ;
  }
  printf( "\n" ) ;
#endif

  gsl_linalg_SV_decomp( alpha , Q , S , Work ) ;

  for( i = 0 ; i < Nstates ; i++ ) {
    //printf( "Singular values %zu %e \n" , i , gsl_vector_get( S , i ) ) ;
    if( gsl_vector_get( S , i ) < gsl_vector_get( S , 0 )*1E-6 ) {
      gsl_vector_set( S , i , 0.0 ) ;
    }
  }
  //exit(1) ;
  
  gsl_linalg_SV_solve( alpha , Q , S , beta , delta ) ;

  //printf( "States \n" ) ;
  double coeffs[ Nstates+1 ] ;

  for( i = 0 ; i < Nstates ; i++ ) {
    coeffs[ i ] = gsl_vector_get( delta , i ) ;
    //printf( "%f\n" , gsl_vector_get( delta , i ) ) ;
  }
  coeffs[Nstates] = 1 ;

  gsl_poly_complex_workspace *w =
    gsl_poly_complex_workspace_alloc( Nstates+1 ) ;
  double *z = malloc( 2 * Nstates * sizeof( double ) ) ;
  
  gsl_poly_complex_solve( coeffs , Nstates+1 , w , z ) ;
  
  //
  //printf( "SOLS\n" ) ;
  for( i = 0 ; i < Nstates ; i++ ) {
    x[ i ] = z[ 2*i ] + I * z[ 2*i+1 ] ;
    //printf( "%f %f \n" , creal( x[i] ) , cimag( x[i] ) ) ;
  }
      
  // free temporary Y-data
  free( Y ) ;
  free( z ) ;

  gsl_poly_complex_workspace_free( w ) ;

  // free the gsl vectors
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( Q ) ;
  
  gsl_vector_free( beta ) ;
  gsl_vector_free( delta ) ;
    
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
  double *e = malloc( NSTATES * sizeof( double ) ) ;
  
  size_t t , i , j ;
  for( t = 0 ; t < smallt ; t++ ) {

    blackbox2( x , data , NDATA , NSTATES , t , NSTATES ) ;
    
    // loggy log log .... x[i] = exp( -e[i] )
    for( i = 0 ; i < NSTATES ; i++ ) {
      e[ i ] = fabs( creal( clog( x[i] ) ) ) ;
      //printf( "%f %f -> %f\n" , creal( x[i] ) , cimag( x[i] ) , e[i] ) ;
    }
    
    // sort the e's ?
    qsort( e , NSTATES , sizeof(double) , comp ) ;
	
    for( i = 0 ; i < NSTATES ; i++ ) {
      masses[i][t] = e[i] ;
    }
  }

  free( e ) ;
  free( x ) ;

  return ;
}
