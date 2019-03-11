/**
   @file inverse.c
   @brief given some distributions it inverts them
 */
#include "gens.h"

#include "cruel_runnings.h"
#include "resampled_ops.h"
#include "stats.h"

#define BD

int
fit_inverse( struct input_params *Input )
{
  size_t i , j ;

  if( Input -> Data.Nsim != 25 ) {
    printf( "Nsim not as expected %zu \n" , Input -> Data.Nsim ) ;
    return FAILURE ;
  }
  
  gsl_matrix *alpha = gsl_matrix_alloc( 5 , 5 ) ;
  gsl_matrix *ainv  = gsl_matrix_alloc( 5 , 5 ) ;
  gsl_permutation *perm = gsl_permutation_alloc( 5 ) ;
  int signum ;
  
  // invert each case by LU
  size_t l , k ;
  printf( "Zinverse matrix\n" ) ;
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {
      printf( "{%e %e} " , Input -> Data.y[ k + 5*l ].avg ,
	      Input -> Data.y[ k + 5*l ].err ) ;
    }
    printf( "\n" ) ;
  }
  printf( "\n" ) ;
  
  for( i = 0 ; i < Input -> Data.y[0].NSAMPLES ; i++ ) {

    #ifdef BD
    for( l = 0 ; l < 5 ; l++ ) {
      for( k = 0 ; k < 5 ; k++ ) {
	gsl_matrix_set( alpha , l , k , 0.0 ) ;
      }
    }
    gsl_matrix_set( alpha , 0 , 0 , Input -> Data.y[0].resampled[i] ) ;
    
    gsl_matrix_set( alpha , 1 , 1 , Input -> Data.y[6].resampled[i] ) ;
    gsl_matrix_set( alpha , 1 , 2 , Input -> Data.y[7].resampled[i] ) ;
    gsl_matrix_set( alpha , 2 , 1 , Input -> Data.y[11].resampled[i] ) ;
    gsl_matrix_set( alpha , 2 , 2 , Input -> Data.y[12].resampled[i] ) ;

    gsl_matrix_set( alpha , 3 , 3 , Input -> Data.y[18].resampled[i] ) ;
    gsl_matrix_set( alpha , 3 , 4 , Input -> Data.y[19].resampled[i] ) ;
    gsl_matrix_set( alpha , 4 , 3 , Input -> Data.y[23].resampled[i] ) ;
    gsl_matrix_set( alpha , 4 , 4 , Input -> Data.y[24].resampled[i] ) ;
    
    #else
    for( l = 0 ; l < 5 ; l++ ) {
      for( k = 0 ; k < 5 ; k++ ) {
	gsl_matrix_set( alpha , l , k ,
			Input -> Data.y[ k + 5*l ].resampled[i] ) ;
      }
    }
    #endif
    // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
    if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
      printf( "[LM] LU decomp broken?\n" ) ;
      return !GSL_SUCCESS ;
    }
    if( gsl_linalg_LU_invert( alpha , perm , ainv ) != GSL_SUCCESS ) {
      printf( "[LM] LU solve broken?\n" ) ;
      return !GSL_SUCCESS ;
    }
    // poke it into the x direction
    for( l = 0 ; l < 5 ; l++ ) {
      for( k = 0 ; k < 5 ; k++ ) {
	Input -> Data.x[ k + 5*l ].resampled[i] =
	  gsl_matrix_get( ainv , l , k ) ;
      }
    }
    //
  }

  printf( "\nZ matrix\n" ) ;
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {
      compute_err( &Input -> Data.x[ k + 5*l ] ) ;
      if( k != 4 ) {
	printf( "%f(%f) & " ,
		Input -> Data.x[ k + 5*l ].avg ,
		Input -> Data.x[ k + 5*l ].err ) ;
      } else {
	printf( "%f(%f) \\\\ \n" ,
		Input -> Data.x[ k + 5*l ].avg ,
		Input -> Data.x[ k + 5*l ].err ) ;
      }
    }
  }

  // write out a flat file
  FILE *file = fopen( "Z.flat" , "w" ) ;
  fprintf( file , "%d\n" , Input -> Data.x[0].restype ) ;
  fprintf( file , "%d\n" , 25 ) ;
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {
      fprintf( file , "%zu\n" , Input -> Data.x[k+5*l].NSAMPLES ) ;
      for( j = 0 ; j < Input -> Data.x[k+5*l].NSAMPLES ; j++ ) {
	fprintf( file , "%1.15e %1.15e\n" , (double)k+5*l , Input -> Data.x[k+5*l].resampled[j] ) ;
      }
    }
  }

  // get alpha_s / ( 4 * M_PI )
  const double a_4pi = run_MZ_2nf3( 0.1182 , Input -> Traj[0].Fit_Low , 5 ) / ( 4*M_PI ) ;

  // convert to MS
  const double dGGMS[5][5] = { { 1. + a_4pi * 0.2118441111462292 , 0 , 0 , 0 , 0 } ,
			       { 0 , 1. + a_4pi *0.04318883230477005 , a_4pi * 0.25913299382862 , 0 , 0 } ,
			       { 0 , a_4pi * 0.8068528194400547 , 1. + a_4pi * 4.495606258202168 , 0 , 0 } ,
			       { 0 , 0 , 0 , 1. + a_4pi * 1.026125409459049 , a_4pi * 0.495776772984037 } ,
			       { 0 , 0 , 0 , - a_4pi * 2.877724730449623 , 1. + a_4pi * 5.737215148874458 } } ;

  printf( "MS-bar conversion matrix\n" ) ;
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {
      printf( "%f " , dGGMS[l][k] ) ;
    }
    printf( "\n" ) ;
  }

  struct resampled *Zms = malloc( 25 * sizeof( struct resampled ) ) ;
  for( l = 0 ; l < 25 ; l++ ) {
    Zms[l] = init_dist( NULL , Input -> Data.x[l].NSAMPLES , Input -> Data.x[l].restype ) ;
  }
  
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {					
      for( j = 0 ; j < 5 ; j++ ) {
	// R[lj] * Z[jk]
	rapby( &Zms[k+5*l] , Input -> Data.x[k+j*5] , dGGMS[l][j] ) ;
      }
      printf( "(%f %f) " , Zms[k+5*l].avg , Zms[k+5*l].err ) ;
    }
    printf( "\n" ) ;
  }

  
  file = fopen( "ZMS_gg.flat" , "w" ) ;
  fprintf( file , "%d\n" , Zms[0].restype ) ;
  fprintf( file , "%d\n" , 25 ) ;
  for( l = 0 ; l < 5 ; l++ ) {
    for( k = 0 ; k < 5 ; k++ ) {
      fprintf( file , "%zu\n" , Zms[k+l*5].NSAMPLES ) ;
      for( j = 0 ; j < Zms[k+5*l].NSAMPLES ; j++ ) {
	fprintf( file , "%1.15e %1.15e\n" , (double)k+5*l , Zms[k+5*l].resampled[j] ) ;
      }
    }
  }

  for( l = 0 ; l < 25 ; l++ ) {
    free( Zms[l].resampled ) ;
  }
  free( Zms ) ;

  fclose( file ) ;
  
  return SUCCESS ;
}
