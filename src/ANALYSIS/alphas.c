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

#include <stdlib.h>

//#define FITRANGES

// qsort comparison function for the bootstrap
int 
comp2( const void *elem1 , 
      const void *elem2 ) 
{

  const double f = *( (double*)elem1 ) ;
  const double s = *( (double*)elem2 ) ;
  if (f > s) { return  1 ; }
  if (f < s) { return -1 ; }
  return 0 ;
}

int
fit_alphas( struct input_params *Input )
{
  test_running( ) ;

  const double mu = 2.00 ;

  if( Input -> Fit.Fitdef == ALPHA_D0_MULTI ) {
    set_mu_multi( mu ) ;
  } else {
    set_mu( mu ) ;
  }

  // compute delta == ( \Pi( q_1^2 ) - \Pi( q^2 ) ) / ( t1 - t2 )
  
  // highest
  //const double q1ref[3] = { 2.0*2.0 , 2.1*2.1 , 2.2*2.2 } ;
  
  // hi
  //const double q1ref[3] = { 1.9*1.9 , 2.0*2.0 , 2.1*2.1 } ;
  //const double q1ref[3] = { 1.8*1.8 , 1.95*1.95 , 2.1*2.1 } ;

  // mid
  const double q1ref[3] = { 1.8*1.8 , 1.9*1.9 , 2.0*2.0 } ;
  
  // low
  //const double q1ref[3] = { 1.7*1.7 , 1.8*1.8 , 1.9*1.9 } ;

  //const double q1ref[3] = { 1.6*1.6 , 1.6*1.6 , 1.6*1.6 } ;


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
    }

    free( tmp.resampled ) ;
    
    shift += Input -> Data.Ndata[i] ;
  }

#ifdef FITRANGES
  // loop lower and upper fit bounds
  double low = 3.7 , fit_super , fit_fine , fit_coarse = 0.0 ;

  size_t Nlow = 0 ;

  struct resampled fullweightave , fullave ;
  fullweightave = init_dist( NULL , Input -> Data.y[0].NSAMPLES ,
			     Input -> Data.y[0].restype ) ;
  fullave = init_dist( NULL , Input -> Data.y[0].NSAMPLES ,
		       Input -> Data.y[0].restype ) ;

  size_t Naccept = 0 ;
  double result[ 1024 ] = {} ;
  
  for( low = 5.5 ; low < 6.75 ; low += 0.25 ) {

    struct resampled aveweight , ave ;
    aveweight = init_dist( NULL , Input -> Data.y[0].NSAMPLES ,
			   Input -> Data.y[0].restype ) ;
    ave = init_dist( NULL , Input -> Data.y[0].NSAMPLES ,
		     Input -> Data.y[0].restype ) ;
    double sumweights = 0.0 ;
    
    //for( fit_coarse = low+0.5 ; fit_coarse < 9 ; fit_coarse+=0.5 ) { 
      for( fit_fine = low+2 ; fit_fine < 11 ; fit_fine+=0.5 ) {
	for( fit_super = fit_fine+2 ; fit_super < 16 ; fit_super+=1.0 ) {

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
	      Input -> Traj[k].Fit_Low = low ;
	      Input -> Traj[k].Fit_High = fit_coarse ;
	      break ;
	    }
	  }

	  double chisq ;
	  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
      
	  struct resampled amz = run_distribution_nf3_2MZ( fit[0] , mu , 4 ) ;

	  if( chisq < 1 ) {
	    printf( "%f %f %f %f ( chi %f ) alpha(%f) -> amz :: %f %f \n" ,
		    low ,
		    fit_super , fit_fine , fit_coarse ,
		    chisq , mu ,
		    amz.avg , amz.err ) ;

	    printf( "CORRECTIONS :: CHI %f | %f +/- %f | %f +/- %f | %f +/- %f \n" ,
		    chisq , 
		    fit[1].avg , fit[1].err ,
		    fit[2].avg , fit[2].err ,
		    fit[3].avg , fit[3].err ) ;

	    // compute alpha
	    add( &ave , amz ) ;

	    result[Naccept] = amz.avg ;

	    // compute weighted alpha
	    mult_constant( &amz , 1.0/chisq ) ;
	    sumweights += 1.0/chisq ;
	
	    add( &aveweight , amz ) ;
	    
	    Naccept++ ;
	  }

	  free_fitparams( fit , Input -> Fit.Nlogic ) ;
	
	  free( amz.resampled ) ;
	}
      }
      //}
      
    divide_constant( &aveweight , sumweights ) ;
    divide_constant( &ave , Naccept ) ;
    
    printf( "------------------------------------\n" ) ;

    printf( "Weighted amz (low,val,err) :: %f %f %f \n" ,
	    low , aveweight.avg , aveweight.err ) ;

    printf( "UnWeighted amz (low,val,err) :: %f %f %f \n" ,
	    low , ave.avg , ave.err ) ;

    add( &fullweightave , aveweight ) ;
    add( &fullave , ave ) ;

    printf( "------------------------------------\n" ) ;

    // free the weighted average
    free( aveweight.resampled ) ;
    free( ave.resampled ) ;

    Nlow++ ;
  }

  divide_constant( &fullweightave , Nlow ) ;
  divide_constant( &fullave , Nlow ) ;

  printf( "FullWeighted (val,err) :: %f %f \n" ,
	  fullweightave.avg , fullweightave.err ) ;

  printf( "FullAve (val,err) :: %f %f \n" ,
	  fullave.avg , fullave.err ) ;

  free( fullweightave.resampled ) ;
  free( fullave.resampled ) ;

  // sort the result
  qsort( result , Naccept , sizeof(double) , comp2 ) ;

  // get the +/- sigma
  const double confidence = 68.2689492 ;
  const double omitted = 0.5 * ( 100. - confidence ) ;
  const size_t bottom = (size_t)( omitted * Naccept / 100. ) ;
  const size_t top = ( Naccept - 1 - bottom ) ;

  printf( "Systematic :: %f %f \n" ,
	  result[top] - 0.117975 ,
	  0.117975 - result[bottom] ) ;

  printf( "Systematic2 :: %f %f %f \n" ,
	  result[Naccept/2] ,
	  result[top] - result[Naccept/2] ,
	  result[Naccept/2] - result[bottom] ) ;

#else

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
      
  struct resampled amz = run_distribution_nf3_2MZ( fit[0] , mu , 4 ) ;

  printf( "%f alpha(%f) -> amz :: %f %f \n" ,
	  chisq , mu , amz.avg , amz.err ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
	
  free( amz.resampled ) ;
  
#endif

  return SUCCESS ;
}
