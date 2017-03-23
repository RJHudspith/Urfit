/**
   @file alphas.c
   @brief driver for computing alphas
 */
#include "gens.h"

#include "alpha_D0.h"
#include "alpha_D0_multi.h"
#include "cruel_runnings.h"
#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"

//#define FITRANGES

int
fit_alphas( struct input_params *Input )
{
  test_running( ) ;

  const double mu = 2.0 ;
  set_mu_multi( mu ) ;
  
  // compute delta == ( \Pi( q_1^2 ) - \Pi( q^2 ) ) / ( t1 - t2 )
  
  // hi2
  //const double q1ref[3] = { 1.75*1.75 , 1.85*1.85 , 1.95*1.95 } ;

  // long
  //const double q1ref[3] = { 1.6*1.6 , 1.75*1.75 , 1.9*1.9 } ;

  // hi
  //const double q1ref[3] = { 1.7*1.7 , 1.8*1.8 , 1.9*1.9 } ;

  // mid
  //const double q1ref[3] = { 1.6*1.6 , 1.7*1.7 , 1.8*1.8 } ;

  const double q1ref[3] = { 1.7*1.7 , 1.7*1.7 , 1.7*1.7 } ;

  // low
  //const double q1ref[3] = { 1.6*1.6 , 1.675*1.675 , 1.75*1.75 } ;
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    double x = 10000 ;
    size_t idx = shift ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      if( fabs( q1ref[i%3] - Input -> Data.x[j].avg ) < x ) {
	x = fabs( q1ref[i%3] - Input -> Data.x[j].avg ) ;
	idx = j ;
      }
    }
    printf( "Q1ref %f (CLOSEST) %f \n" , q1ref[i%3] , Input -> Data.x[idx].avg ) ;

    if( Input -> Fit.Fitdef == ALPHA_D0_MULTI ) {
      set_Q1_multi( Input -> Data.x[idx].avg , i ) ;
    } else if( Input -> Fit.Fitdef == ALPHA_D0 ) {
      set_Q1( Input -> Data.x[idx].avg , i ) ;
    }
    
    // now we do the subtraction
    struct resampled tmp = init_dist( &Input -> Data.y[idx] ,
				      Input -> Data.y[idx].NSAMPLES ,
				      Input -> Data.y[idx].restype ) ;

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      
      subtract( &Input -> Data.y[j] , tmp ) ;

      if( j != idx ) {
	divide_constant( &Input -> Data.y[j] ,
			 log( Input -> Data.x[idx].avg / ( Input -> Data.x[j].avg ) ) ) ;
      }
      
      mult_constant( &Input -> Data.y[j] , -4 * M_PI * M_PI ) ;
      
      subtract_constant( &Input -> Data.y[j] , 1.0 ) ;

      printf( "%f %f %f\n" , Input -> Data.x[j].avg , Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
    }
    printf( "\n" ) ;
    
    free( tmp.resampled ) ;
    
    shift += Input -> Data.Ndata[i] ;
  }
  //exit(1) ;

#ifdef FITRANGES
  // loop lower and upper fit bounds
  double low = 3.7 , fit_super , fit_fine ;
  
  for( low = 4.0 ; low < 5.2 ; low += 0.2 ) {

    struct resampled aveweight ;
    aveweight = init_dist( NULL , Input -> Data.y[0].NSAMPLES ,
			   Input -> Data.y[0].restype ) ;
    double sumweights = 0.0 ;
    
    for( fit_fine = 6 ; fit_fine < 13 ; fit_fine++ ) {
      for( fit_super = fit_fine+1 ; fit_super < 22 ; fit_super++ ) {

	size_t k ;
	for( k = 0 ; k < 12*3 ; k++ ) {
	  switch( k/12 ) {
	  case 0 :
	    Input -> Traj[k].Fit_Low = low ;
	    Input -> Traj[k].Fit_High = fit_super ;
	    break ;
	  case 1 :
	    Input -> Traj[k].Fit_Low = low ;
	    Input -> Traj[k].Fit_High = fit_fine ;
	    break ;
	  case 2 :
	    break ;
	  }
	}

	double chisq ;
	struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
      
	struct resampled amz = run_distribution_nf3_2MZ( fit[0] , mu ) ;

	printf( "%f %f %f ( chi %f ) alpha(%f) -> amz :: %f %f \n" ,
		low , fit_super , fit_fine ,
		chisq , mu ,
		amz.avg , amz.err ) ;

	printf( "CORRECTIONS :: CHI %f | %f +/- %f | %f +/- %f | %f +/- %f | CHI :: %f \n" ,
		chisq , 
		fit[1].avg , fit[1].err ,
		fit[2].avg , fit[2].err ,
		fit[3].avg , fit[3].err ) ;

	// compute alpha
	mult_constant( &amz , 1.0/chisq ) ;
	sumweights += 1.0/chisq ;
	
	add( &aveweight , amz ) ;

	free_fitparams( fit , Input -> Fit.Nlogic ) ;
	
	free( amz.resampled ) ;
      }
    }
    divide_constant( &aveweight , sumweights ) ;

    printf( "------------------------------------\n" ) ;

    printf( "Weighted amz (low,val,err) :: %f %f %f \n" ,
	    low , aveweight.avg , aveweight.err ) ;

    printf( "------------------------------------\n" ) ;

    // free the weighted average
    free( aveweight.resampled ) ;
  }
#else

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
      
  struct resampled amz = run_distribution_nf3_2MZ( fit[0] , mu ) ;

  printf( "%f alpha(%f) -> amz :: %f %f \n" ,
	  chisq , mu , amz.avg , amz.err ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
	
  free( amz.resampled ) ;
#endif
  
  return SUCCESS ;
}
