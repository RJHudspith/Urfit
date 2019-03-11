/**
   @file adler.c
   @brief adler function analysis
 */
#include "gens.h"

#include "adler_alpha_D0.h"
#include "adler_alpha_D0_multi.h"
#include "cruel_runnings.h"
#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"
#include "rng.h"
#include "stats.h"
#include "write_flat.h"

// perform an SVD to get the continuum result
//
// Matrix is of the form :
// | 1 bar{mpi^2} bar{mk^2} a^2 | | y_0 |   | D(0) |
// | 1 bar{mpi^2} bar{mk^2} a^2 | | C_1 | = | D(1) |
// | 1 bar{mpi^2} bar{mk^2} a^2 | | C_2 | = | D(2) |
// | 1 bar{mpi^2} bar{mk^2} a^2 | | C_a | = | D(3) |
// |  ... .... ... ... ... ...  |           | D(.) |
static int
get_cont( struct input_params *Input )
{
  // construct the matrix M columns N rows
  const size_t M = 2 , N = 2 ;
  
  if( Input -> Data.Nsim < N ) {
    printf( "Not enough data points %zu, exiting\n" , Input -> Data.Nsim ) ;
    return FAILURE ;
  }
  
  // and finally the lattice spacings squared
  const double a2[ 3 ] = { 0.04949440885942718 ,
			   0.07684254741528086 ,
			   0.16601189178800163 } ;

  // construct the matrix M columns M rows as a GLS fit
  gsl_matrix *alpha = gsl_matrix_alloc( M , M ) ;
  gsl_matrix *Q     = gsl_matrix_alloc( M , M ) ;
  gsl_vector *S     = gsl_vector_alloc( M ) ;
  gsl_vector *Work  = gsl_vector_alloc( M ) ;
  
  gsl_vector *beta  = gsl_vector_alloc( M ) ;
  gsl_vector *delta = gsl_vector_alloc( M ) ;

  double **U = malloc( N * sizeof( double ) ) ;

  // do the matrix stuff, only need to invert alpha once!
  size_t i ;
  for( i = 0 ; i < N ; i++ ) {
    U[ i ] = malloc( M * sizeof( double ) ) ;
    U[i][0] = 1 ;
    U[i][1] = a2[i] ;
  }

  struct resampled *c0 = malloc( Input -> Data.Ndata[0] * sizeof( struct resampled ) ) ;
  struct resampled *c1 = malloc( Input -> Data.Ndata[0] * sizeof( struct resampled ) ) ;
  double *Weight = calloc( N , sizeof( double ) ) ;

  // compute the correct normalisation
  double NORM = 1.0 ;
  const size_t NSAMPLES = Input -> Data.x[0].NSAMPLES ;

  // get the norm correct
  switch( Input -> Data.x[0].restype ) {
  case Raw :
    NORM = 1.0 / (double)( NSAMPLES * ( NSAMPLES - 1.0 ) ) ;
    break ;
  case JackKnife :
    NORM = ( 1.0 - 1.0/(double)NSAMPLES ) ;
    break ;
  case BootStrap :
    NORM = 1.0 / (double)NSAMPLES ;
    break ;
  }
  
  size_t j , k ;
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    
    c0[j] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;
    c1[j] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;
    
    size_t i , l , m ;
    for( i = 0 ; i < N ; i++ ) {    
      const register double ave = Input -> Data.y[ j + i*Input -> Data.Ndata[0] ].avg ;
      register double sum = 0.0 ;
      size_t k ;
      for( k = 0 ; k < NSAMPLES ; k++ ) {
	sum += 
	  ( Input -> Data.y[ j + i*Input -> Data.Ndata[0] ].resampled[k] - ave ) *  
	  ( Input -> Data.y[ j + i*Input -> Data.Ndata[0] ].resampled[k] - ave ) ;
      }
      Weight[i] = 1.0 / ( sum * NORM ) ;
    }

    // set the alpha matrix
    for( i = 0 ; i < M ; i++ ) {
      for( l = 0 ; l < M ; l++ ) {
	double sum = 0.0 ;
	for( m = 0 ; m < N ; m++ ) {
	  sum += U[m][i] * Weight[m] * U[m][l] ;
	}
	gsl_matrix_set( alpha , i , l , sum ) ;
      }
    }
    
    // svds
    if( gsl_linalg_SV_decomp( alpha , Q , S , Work ) != GSL_SUCCESS ) {
      printf( "SVD comp failure \n" ) ;
      return FAILURE ;
    } 
		
    for( k = 0 ; k < Input -> Data.x[0].NSAMPLES ; k++ ) {

      // set the solution vector U^T y
      for( i = 0 ; i < M ; i++ ) {
	double sum = 0.0 ;
	for( m = 0 ; m < N ; m++ ) {
	  sum += U[m][i] * Weight[m] * Input -> Data.y[ j + m*Input -> Data.Ndata[0] ].resampled[ k ] ;
	}
	gsl_vector_set( beta , i , sum ) ;
      }
      
      // do a solve
      if( gsl_linalg_SV_solve( alpha , Q , S , beta , delta ) != GSL_SUCCESS ) {
	return FAILURE ;
      }
      // poke delta into our struct
      c0[j].resampled[k] = gsl_vector_get( delta , 0 ) ;
      c1[j].resampled[k] = gsl_vector_get( delta , 1 ) ;
    }
    // do the average
    for( i = 0 ; i < M ; i++ ) {
      double sum = 0.0 ;
      for( m = 0 ; m < N ; m++ ) {
	sum += U[m][i] * Weight[m] * Input -> Data.y[ j + m*Input -> Data.Ndata[i] ].avg ;
      }
      gsl_vector_set( beta , i , sum ) ;
    }
    if( gsl_linalg_SV_solve( alpha , Q , S , beta , delta ) != GSL_SUCCESS ) {
      return FAILURE ;
    }
    // poke delta into our struct
    c0[j].avg = gsl_vector_get( delta , 0 ) ;
    c1[j].avg = gsl_vector_get( delta , 1 ) ;
    
    compute_err( &c0[j] ) ;
    compute_err( &c1[j] ) ;
    
    printf( "ADLERCONT %e %e %e\n" , Input -> Data.x[j].avg , c0[j].avg , c0[j].err ) ;
    printf( "FRAC %e \n" , 100*c0[j].err/c0[j].avg ) ;
    printf( "CA %e %e %e\n" , Input -> Data.x[j].avg , c1[j].avg , c1[j].err ) ;
    
    // compute a chi^2
    double chi = 0.0 ;
    for( i = 0 ; i < N ; i++ ) {
      double loc_sum = 0.0 ;
      for( m = 0 ; m < M ; m++ ) {
	loc_sum += U[i][m] * gsl_vector_get( delta , m ) ;
      }
      chi +=
	( Input -> Data.y[ j + i*Input -> Data.Ndata[0] ].avg - loc_sum ) *
	Weight[ i ] *
	( Input -> Data.y[ j + i*Input -> Data.Ndata[0] ].avg - loc_sum ) ;
    }
    printf( "CHI/DOF :: %e %e \n" , Input -> Data.x[j].avg , chi / N ) ;

    // overwrite Y?
  } 

  for( i = 0 ; i < Input -> Data.Ndata[0] ; i++ ) {
    free( c0[i].resampled ) ;
    free( c1[i].resampled ) ;
  }
  free( c0 ) ;
  free( c1 ) ;
  
  free( Weight ) ;

  for( i = 0 ; i < N ; i++ ) {
    free( U[i] ) ;
  }
  free( U ) ;

  // free the gsl vectors
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( Work ) ;
  
  gsl_vector_free( beta ) ;
  gsl_vector_free( delta ) ;

  return 0 ;
}

const static double
lp1( const double t , const double d5 ) {
  return 1.0 ;
}

const static double
lp2( const double t , const double d5 ) {
  return 1.63982 - 2.25 * t ;
}

const static double
lp3( const double t , const double d5 ) {
  return 6.37108 - 11.379195 * t + 5.0625 * t ;
}

const static double
lp4( const double t , const double d5 ) {
  return 49.0769 -66.182813*t + 47.404784*t*t -11.390625*t*t*t ;
}

const static double
lp5( const double t , const double d5 ) {
  return d5 -598.354375*t + 388.732597*t*t - 162.464353*t*t*t + 25.628906*t*t*t*t ;
}

const static double (*loop[5])( const double t , const double d5 ) = { lp1 , lp2 , lp3 , lp4 , lp5 } ;

// perturbation theory result
double
PT( const double alpha ,
    const double mu ,
    const double Q )
{
  const double t = log( Q*Q / ( mu*mu ) ) ;
  const double a_pi = alpha/M_PI ;
  register double pt = 0.0 ;
  size_t loops ;
  for( loops = 5 ; loops > 0 ; loops-- ) {
    pt = a_pi * ( loop[ loops-1 ](t,400) + pt ) ;
  }
  return pt ;
}


#define SINGLE

int
adler_analysis( struct input_params *Input )
{
#ifdef MULTI
  if( Input -> Data.Nsim > 12 ) {
    fprintf( stderr , "New analysis only expects 12 averaged bootstraps!\n" ) ;
    exit(1) ;
  }

  const double mu = 3.0 ;

  set_mu_multi_adler( mu ) ;

  /*
  const double ainv[ 9 ] = { 4.494919612756264 , 4.494919612756264 , 4.494919612756264 ,
			     3.607440054844607 , 3.607440054844607 , 3.607440054844607 ,
			     2.454315559701493 , 2.454315559701493 , 2.454315559701493 } ;
  */

  const size_t ZMAP[ 9 ] = { 0 , 0 , 0 , 
			     1 , 1 , 1 ,
			     2 , 2 , 2 } ;

#else
  
  const double mu = 2.00 ;

  set_mu_adleralpha( mu ) ;
  
  //const double ainv[ 3 ] = { 4.494919612756264 , 3.607440054844607 , 2.454315559701493 } ;
  //const double ainv[ 3 ] = { 2.454315559701493 } ;
  const size_t ZMAP[ 3 ] = { 0 , 1 , 2 } ;
#endif
  
  // create a fake gaussian distribution
  struct resampled *Z = malloc( 3 * sizeof( struct resampled ) ) ;

  Z[0] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;
  Z[1] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;
  Z[2] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;

  size_t k ;
  Z[0].avg = 0.9699 ;
  Z[1].avg = 0.9636 ;
  Z[2].avg = 0.9553 ;

  printf( "%e %e \n" , Z[0].avg , Z[0].err ) ;
  printf( "%e %e \n" , Z[1].avg , Z[1].err ) ;
  printf( "%e %e \n" , Z[2].avg , Z[2].err ) ;

  init_rng( 1234567 ) ;
  
  for( k = 0 ; k < Z[0].NSAMPLES ; k++ ) {
    Z[0].resampled[k] = 0.9699 ;//+ rng_gaussian( 0.0026 ) ;
  }
  for( k = 0 ; k < Z[0].NSAMPLES ; k++ ) {
    Z[1].resampled[k] = 0.9636 ;//+ rng_gaussian( 0.0034 ) ;
  }
  for( k = 0 ; k < Z[0].NSAMPLES ; k++ ) {
    Z[2].resampled[k] = 0.9553 ;//+ rng_gaussian( 0.0053 ) ;
  }
  raise( &Z[0] , 2 ) ;
  raise( &Z[1] , 2 ) ;
  raise( &Z[2] , 2 ) ;

  free_rng() ;

  printf( "%e %e \n" , Z[0].avg , Z[0].err ) ;
  printf( "%e %e \n" , Z[1].avg , Z[1].err ) ;
  printf( "%e %e \n" , Z[2].avg , Z[2].err ) ;
     
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      
      // mult by Z_V^2
      mult( &( Input -> Data.y[j] ) , Z[ ZMAP[i] ] ) ;
      
      // multiply by 4\pi^2
      mult_constant( &( Input -> Data.y[j] ) , (4.*M_PI*M_PI) ) ;

      // subtract 1
      subtract_constant( &( Input -> Data.y[j] ) , 1.0 ) ;

      #if 0
      size_t k ;
      char str[256] ;
      sprintf( str , "./BOOTS/bootdist_%zu", j ) ;
      FILE *dist = fopen( str , "w" ) ;
      for( k = 0 ; k < Input -> Data.y[j].NSAMPLES ; k++ ) {
	fprintf( dist , "%e\n" , Input -> Data.y[j].resampled[k] ) ;
      }
      fclose( dist ) ;
      #endif
      
      // divide by Q^2
      //divide( &( Input -> Data.y[j] ) , Input -> Data.x[j] ) ;

      printf( "%e %e %e \n" , Input -> Data.x[j].avg , Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
    }
    printf( "\n" ) ;
    shift += Input -> Data.Ndata[i] ;
  }

  free( Z[0].resampled ) ;
  free( Z[1].resampled ) ;
  free( Z[2].resampled ) ;
  free( Z ) ;

  get_cont( Input ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  struct resampled amz = run_distribution_nf3_2MZ( Fit[0] , mu , 4 ) ;

  printf( "%f alpha(%f) -> amz :: %f %f \n" ,
	  Chi , mu , amz.avg , amz.err ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
	
  free( amz.resampled ) ;
  
  return SUCCESS ;
}


