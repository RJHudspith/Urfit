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

//#define FINE_AVG

#define MID_MS

int
old_adler_analysis( struct input_params *Input )
{
#ifdef MULTI
  if( Input -> Data.Nsim > 12 ) {
    fprintf( stderr , "New analysis only expects 12 averaged bootstraps!\n" ) ;
    exit(1) ;
  }

  const double mu = 2.5 ;

  set_mu_multi_adler( mu ) ;
  
  const size_t ZMAP[ 9 ] = { 0 , 0 , 0 , 
			     1 , 1 , 1 ,
			     2 , 2 , 2 } ;
  
#else
  
  const double mu = 3.00 ;

  set_mu_adleralpha( mu ) ;
   

#endif

  #if (defined COARSE_AVG)
  
  const double ainv[ 4 ] = { 2.454315559701493 , 2.454315559701493 , 2.454315559701493 , 2.454315559701493  } ;
  const double mcorr[ 4 ] = { 1.0141483841066723 , 1.0141483841066723 , 1.0244390171941844 , 1.0244390171941844 } ;
  const double z[2] = { 0.96928 , 0.96947 } ;
  const double dZ[2] = { 0.00028 , 0.00021 } ;

  const size_t ZMAP[ 4 ] = { 0 , 0 , 1, 1 } ;

  #elif (defined MID_AVG)

  const double ainv[ 4 ] = { 3.607440054844607 , 3.607440054844607 , 3.607440054844607 , 3.607440054844607 } ;
  
  const double mcorr[ 4 ] = { 1.008453217915728 , 1.008453217915728 , 1.01619406867846 , 1.01619406867846  } ;
  const double z[2] = { 0.96553 , 0.96596 } ;
  const double dZ[2] = { 0.00012 , 0.00006 } ;

  const size_t ZMAP[ 4 ] = { 0 , 0 , 1, 1 } ;

#elif (defined FINE_AVG)

  const double ainv[ 4 ] = { 4.494919612756264 , 4.494919612756264 , 4.494919612756264 , 4.494919612756264 } ;
  
  const double mcorr[ 4 ] = { 1.0060271084064631 , 1.0060271084064631 , 1.0060271084064631 , 1.0060271084064631 } ;
  
  const double z[2] = { 0.96650 , 0.96650 } ;
  const double dZ[2] = { 0.00011 , 0.00011 } ;
  const size_t ZMAP[ 4 ] = { 0 , 0 , 1, 1 } ;

#elif (defined COARSE_MS)
  
  const double ainv[ 4 ] = { 2.454315559701493 , 2.454315559701493 , 2.454315559701493 , 2.454315559701493  } ;
  const double mcorr[ 4 ] = { 1.0141483841066723 , 1.0141483841066723 , 1.0244390171941844 , 1.0244390171941844 } ;
  const double z[2] = { 0.96928 , 0.96947 } ;
  const double dZ[2] = { 0.00028 , 0.00021 } ;
  const size_t ZMAP[ 4 ] = { 0 , 0 , 1, 1 } ;

#elif (defined MID_MS)

  const double ainv[ 4 ] = { 3.607440054844607 , 3.607440054844607 , 3.607440054844607 , 3.607440054844607 } ;
  const double mcorr[ 4 ] = { 1.008453217915728 , 1.008453217915728 , 1.01619406867846 , 1.01619406867846  } ;
  const double z[2] = { 0.96553 , 0.96596 } ;
  const double dZ[2] = { 0.00012 , 0.00006 } ;
  const size_t ZMAP[ 4 ] = { 0 , 0 , 1, 1 } ;
  
#else

  const double ainv[ 3 ] = { 4.494919612756264 , 3.607440054844607 , 2.454315559701493 } ;
  const double mcorr[ 3 ] = { 1.0060271084064631 , 1.008453217915728 , 1.0141483841066723 } ;
  const double z[3] = { 0.96650 , 0.96553 , 0.96928 } ;
  const double dZ[3] = { 0.00011 , 0.00006 , 0.00028 } ;

  //const size_t ZMAP[ 3 ] = { 0 , 1 , 2 } ;
  #endif

  
  // create a fake gaussian distribution
  struct resampled *Z = malloc( 3 * sizeof( struct resampled ) ) ;

  Z[0] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;
  Z[1] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;
  Z[2] = init_dist( NULL , Input -> Data.x[0].NSAMPLES , Input -> Data.x[0].restype ) ;

  size_t k ;
  Z[0].avg = z[0] ;
  Z[1].avg = z[1] ;
  Z[2].avg = z[2] ;

  printf( "%e %e \n" , Z[0].avg , Z[0].err ) ;
  printf( "%e %e \n" , Z[1].avg , Z[1].err ) ;
  printf( "%e %e \n" , Z[2].avg , Z[2].err ) ;

  init_rng( 1234567 ) ;
  
  for( k = 0 ; k < Z[0].NSAMPLES ; k++ ) {
    Z[0].resampled[k] = Z[0].avg + rng_gaussian( dZ[0] ) ;
  }
  for( k = 0 ; k < Z[0].NSAMPLES ; k++ ) {
    Z[1].resampled[k] = Z[1].avg + rng_gaussian( dZ[1] ) ;
  }
  for( k = 0 ; k < Z[0].NSAMPLES ; k++ ) {
    Z[2].resampled[k] = Z[2].avg + rng_gaussian( dZ[2] ) ;
  }
  raise( &Z[0] , 2 ) ;
  raise( &Z[1] , 2 ) ;
  raise( &Z[2] , 2 ) ;

  free_rng() ;

  printf( "%e %e \n" , Z[0].avg , Z[0].err ) ;
  printf( "%e %e \n" , Z[1].avg , Z[1].err ) ;
  printf( "%e %e \n" , Z[2].avg , Z[2].err ) ;

  struct resampled *y = malloc( Input -> Data.Ndata[0] * sizeof( struct resampled ) ) ;
  struct resampled *x = malloc( Input -> Data.Ndata[0] * sizeof( struct resampled ) ) ;

  for( size_t i = 0 ; i < Input -> Data.Ndata[0] ; i++ ) {
    x[i] = init_dist( NULL , Input -> Data.y[i].NSAMPLES , Input -> Data.y[i].restype ) ;
    y[i] = init_dist( NULL , Input -> Data.y[i].NSAMPLES , Input -> Data.y[i].restype ) ;
  }
     
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    size_t idx=0 ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      mult_constant( &(Input -> Data.y[j]) , mcorr[ i ] ) ;
      
      // multiply x by a^2
      //root( &( Input -> Data.x[j] ) ) ;
      mult_constant( &( Input -> Data.x[j] ) , ainv[i]*ainv[i] ) ;
      
      // mult by Z_V^2
      mult( &( Input -> Data.y[j] ) , Z[ ZMAP[i] ] ) ;
      
      // multiply by 4\pi^2
      mult_constant( &( Input -> Data.y[j] ) , (4.*M_PI*M_PI) ) ;

      // subtract 1
      subtract_constant( &( Input -> Data.y[j] ) , 1.0 ) ;

      add( &y[idx] , Input -> Data.y[j] ) ;
      add( &x[idx] , Input -> Data.x[j] ) ;
      idx++ ;
      
      printf( "%e %e %e \n" , Input -> Data.x[j].avg , Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
    }
    printf( "\n" ) ;
    shift += Input -> Data.Ndata[i] ;
  }

#ifndef COARSE_MS
  FILE *file = fopen( "Fine.flat" , "w" ) ;
  fprintf( file , "2\n" ) ;
  fprintf( file , "%zu\n" , Input -> Data.Ndata[0] ) ;
  for( size_t i = 0 ; i < Input -> Data.Ndata[0] ; i++ ) {
    fprintf( file , "%zu\n" , y[i].NSAMPLES ) ;
    divide_constant( &y[i] , 1.0 ) ;
    divide_constant( &x[i] , 1.0 ) ;
    for( size_t k = 0 ; k < y[i].NSAMPLES ; k++ ) {
      fprintf( file , "%1.15e %1.15e\n" , x[i].resampled[k] , y[i].resampled[k] ) ;
    }	
    fprintf( file , "AVG %1.15e %1.15e\n" , x[i].avg , y[i].avg ) ;
  }
  fclose(file) ;
#endif  

  free( Z[0].resampled ) ;
  free( Z[1].resampled ) ;
  free( Z[2].resampled ) ;
  free( Z ) ;

    // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  /*
  struct resampled amz = run_distribution_nf3_2MZ( Fit[0] , mu , 4 ) ;
  printf( "%f alpha(%f) -> amz :: %f %f \n" ,
	  Chi , mu , amz.avg , amz.err ) ;
  free( amz.resampled ) ;
  */
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
}

int
adler_analysis( struct input_params *Input )
{
  const double mu = 2.6 ;

  set_mu_multi_adler( mu ) ;
  set_mu_adleralpha( mu ) ;
  
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


