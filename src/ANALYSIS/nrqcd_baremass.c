/**
   @file nrqcd_baremass.c
   @brief computes the bare mass from a straight line
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"

int
nrqcd_baremass_analysis( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  const double ainv = 2.206 ; //2.194 ;

  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      printf( "MKIN %f %f %f\n" ,
	      Input -> Data.x[j].avg ,
	      Input -> Data.y[j].avg ,
	      Input -> Data.y[j].err ) ;
      mult_constant( &Input -> Data.y[j] , ainv ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }

  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    spin_average( &Input -> Data.y[j] , Input -> Data.y[j+Input->Data.Ndata[0]] ) ;
    //divide_constant( &Input -> Data.y[j] , 9.445 ) ;
  }

  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  const double Meta = 9.445 ;

  struct resampled res = init_dist( NULL , Input -> Data.y[0].NSAMPLES ,
				    Input -> Data.y[0].restype ) ;
  struct resampled savg = init_dist( &fit[1] ,
				     Input -> Data.y[0].NSAMPLES ,
				     Input -> Data.y[0].restype ) ;

  equate_constant( &res , Meta ,
		   Input -> Data.y[0].NSAMPLES ,
		   Input -> Data.y[0].restype ) ;

  subtract( &res , fit[0] ) ;
  divide( &res , fit[1] ) ;
  
  fprintf( stdout , "Pred %f +/- %f \n" , res.avg , res.err ) ; 

  free( res.resampled ) ;


  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
