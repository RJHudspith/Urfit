/**
   @file autocorr.c
   @brief autocorrelation

   This code computes the integtated autocorrelation time :

   It expects a file of length N double precision measurements

   Defining c(t) = ( Y(t) - \bar{Y} )

   C(T) = \sum_{t=0}^{N} ( c(t) c(t+T) )

   R(T) = C(T) / C(0)

   Setting a cutoff point "n"

   Tau(n) = 0.5 + Nsep * \sum_{T}^{n} R(T)

   We estimate the error on Tau(T) by

   S_E(n) = n * \sqrt( ( 0.5 + \sum_{t}^{n} R(T) ) / N ) 

   The computation of C(T) is performed by convolution with FFTs

   The autocorrelation time is written to the file specified
 */
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "gens.h"
#include "stats.h"

#define HAVE_FFTW3_H

#ifdef HAVE_FFTW3_H

#include <fftw3.h>

// some defines
#define message(a)( fprintf( stdout , "--> %s <--\n\n" , a ) )

// if we have these 
#if ( defined OMP_FFTW ) && ( defined HAVE_OMP_H )
 #include <omp.h>
 static int nthreads = 1 ;
#endif

// for ease of reading
enum{ NOPLAN = 0 } ;

// have a look for parallel ffts
static int
parallel_ffts( void )
{
  // initialise parallel fft ?
#ifdef OMP_FFTW
  if( fftw_init_threads( ) == 0 ) {
    return FAILURE ;
  } else {
    // in here I set the number of fftw threads to be the same
    // as the usual number of parallel threads ..
    #pragma omp parallel
    { nthreads = omp_get_num_threads( ) ; } // set nthreads
  }
  fftw_plan_with_nthreads( nthreads ) ;
  fprintf( stdout , "[PAR] FFTW using %d thread(s) \n" , nthreads ) ;
#endif
  return SUCCESS ;
}

// return 
int
autocorrelation( const struct resampled RAW ,
		 const size_t NSEP ,
		 const char *output )
{
  printf( "In this guy !!! ") ;
  
  // openmp'd fftws
  parallel_ffts( ) ;

  if( RAW.restype != Raw ) {
    fprintf( stderr , "Resampled data is not RAW ..."
	     " Cannot compute autocorrelation\n" ) ;
    return FAILURE ;
  }

  // some constants
  const size_t N = RAW.NSAMPLES ;
  const size_t N2 = 2 * N ;

  fprintf( stdout , "RAWDATA has %zu samples\n" , N ) ;

  fprintf( stdout , "Measurement separation %zu\n\n" , NSEP ) ;

  // allocate memory
  double complex *in  = calloc( N2 , sizeof( double complex ) ) ;
  double complex *out = calloc( N2 , sizeof( double complex ) ) ;

  // subtract the average from each data point
  size_t i ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < N ; i++ ) {
    in[ i ] = ( RAW.resampled[ i ] - RAW.avg ) ;
  }

  // are we doing this using openmp ffts?
#if ( defined OMP_FFTW ) && ( defined HAVE_OMP_H )
  if( parallel_ffts( ) == FAILURE ) {
    printf( "Parallel FFT setting failed \n" ) ;
    return FAILURE ;
  }
#endif
  
  // create the plans
  const fftw_plan forward = fftw_plan_dft_1d( N2 , in , out , 
					      FFTW_FORWARD , FFTW_ESTIMATE ) ; 

  const fftw_plan backward = fftw_plan_dft_1d( N2 , out , in , 
					       FFTW_BACKWARD , FFTW_ESTIMATE ) ;
  
  fftw_execute( forward ) ;

  // convolve
#pragma omp parallel for private(i)
  for( i = 0 ; i < N2 ; i++ ) {
    out[i] = creal( out[i] ) * creal( out[i] ) + 
             cimag( out[i] ) * cimag( out[i] ) ;
  }

  fftw_execute( backward ) ;

  // normalise
  const double zeropoint = 1.0 / in[ 0 ] ;
#pragma omp parallel for private(i)
  for( i = 0 ; i < N2 ; i++ ) {
    in[ i ] *= zeropoint ;
  }

  // summy the lags
  message( "Computing tau(n)" ) ;

  FILE *output_file = fopen( output , "w" ) ;

  printf( "Writing tau(n) to file %s \n" , output ) ;

  size_t n , j ;
  for( n = 0 ; n < N ; n++ ) {
    register double sum = 0.5 ;
    for( j = 0 ; j < n ; j++ ) {
      sum += NSEP * in[j] ;
    }
    // simple error estimate
    const double err = n * sqrt( ( sum ) / N ) ;
    if( isnan( sum ) || isnan( err ) ) break ;
    fprintf( output_file , "%zu %e %e \n" , n * NSEP , sum , err ) ;
  }

  fclose( output_file ) ;

  // memory free
  free( in ) ;
  free( out ) ;
  fftw_destroy_plan(forward ) ;
  fftw_destroy_plan( backward ) ;
#ifdef OMP_FFTW
  // parallel
  fftw_cleanup_threads( ) ;
#endif  
  fftw_cleanup( ) ;

  return SUCCESS ;
}

#else

int
autocorrelation( const struct resampled RAW ,
		 const size_t NSEP ,
		 const char *output )
{
  return SUCCESS ;
}

#endif

// wrapper for the autocorr measurement
int
ACmeasure( const struct input_params Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input.Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input.Data.Ndata[i] ; j++ ) {
      char str[ 64 ] = { 0 } ;
      sprintf( str , "TAUINT_%zu" , j ) ;
      compute_err( &Input.Data.y[j] ) ;
      if( autocorrelation( Input.Data.y[j] ,
			   Input.Traj[i].Increment ,
			   str ) == FAILURE ) {
	return FAILURE ;
      }
    }
    shift = j ;
  }
  return SUCCESS ;
}
