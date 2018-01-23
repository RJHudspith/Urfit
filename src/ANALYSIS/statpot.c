/**
   @file statpot.c
   @brief static potential calculator
 */
#include "gens.h"

#include "fit_and_plot.h"

int
statpot_analysis( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      // sqrt the x
      root( &( Input -> Data.x[j] ) ) ;
      // mult V(r) by x so that we can fit a polynomial
      //mult( &( Input -> Data.y[j] ) , Input -> Data.x[j] ) ;
      //mult( &( Input -> Data.y[j] ) , Input -> Data.x[j] ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  root( &Fit[0] ) ;

  printf( "a\sqrt{sigma} = %f %f \n" , Fit[0].avg , Fit[0].err ) ;

  // write out a flat file
  FILE *file = fopen( "sigma.flat" , "w" ) ;
  fprintf( file , "%zu\n" , Fit[0].NSAMPLES ) ;
  for( j = 0 ; j < Fit[0].NSAMPLES ; j++ ) {
    fprintf( file , "%1.15e %1.15e\n" ,
	     1./(double)Input->Traj[0].Dimensions[3] ,
	     Fit[0].resampled[j] ) ;
  }
  fprintf( file , "AVG %1.15e %1.15e\n" ,
	   1./(double)Input->Traj[0].Dimensions[3] ,
	   Fit[0].avg ) ;
  fclose( file ) ;
  
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}
