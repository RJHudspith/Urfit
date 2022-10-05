/**
   @file nrqcd_baremass.c
   @brief computes the bare mass from a straight line
 */
#include "gens.h"

#include "init.h"
#include "fit_and_plot.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"

int
nrqcd_baremass_analysis( struct input_params *Input )
{
  size_t i , j , shift = 0 ;

#if 0
  const double ainv = 2.194 ; //2.287 ; //2.206 ; //2.194 ;

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
#endif

  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  const double Meta = 4.746473196639064 ;

  struct resampled res = init_dist( NULL , Input -> Data.y[0].NSAMPLES ,
				    Input -> Data.y[0].restype ) ;

  equate_constant( &res , Meta ,
		   Input -> Data.y[0].NSAMPLES ,
		   Input -> Data.y[0].restype ) ;

  subtract( &res , fit[0] ) ;
  divide( &res , fit[1] ) ;
  
  fprintf( stdout , "Pred %1.10f +/- %1.10f \n" , res.avg , res.err ) ;

  raise( &res , -1 ) ;

  fprintf( stdout , "Kc %1.10f +/- %1.10f \n" , res.avg , res.err ) ;

  free( res.resampled ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  /*
  const double ainv = 0.1973269804 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      divide_constant( &Input -> Data.y[j] , ainv ) ;
      raise( &Input -> Data.y[j] , -1 ) ;
      fprintf( stdout , "%f %f\n" , Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  */
  
  return SUCCESS ;
}
