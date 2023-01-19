/**
   @file statpot.c
   @brief static potential calculator
 */
#include "gens.h"

#include "init.h"
#include "fit_and_plot.h"
#include "resampled_ops.h"

//#define MUL_R
//#define MUL_R2

int
statpot_analysis( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      // sqrt the x
      root( &( Input -> Data.x[j] ) ) ;
      // mult V(r) by x so that we can fit a polynomial
      #ifdef MUL_R
      mult( &( Input -> Data.y[j] ) , Input -> Data.x[j] ) ;
      #elif (defined MUL_R2)
      mult( &( Input -> Data.y[j] ) , Input -> Data.x[j] ) ;
      mult( &( Input -> Data.y[j] ) , Input -> Data.x[j] ) ;
      #endif
    }
    shift += Input -> Data.Ndata[i] ;
  }

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  root( &Fit[0] ) ;

  fprintf( stdout , "a sqrt{sigma} = %f %f \n" , Fit[0].avg , Fit[0].err ) ;

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


int
statpot_analysis_v2( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  struct resampled temp = init_dist( NULL ,
				     Input->Data.y[0].NSAMPLES ,
				     Input->Data.y[0].restype ) ;

  FILE *effmass = fopen( "effmass.dat" , "w" ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      // sqrt the x
      root( &( Input -> Data.x[j] ) ) ;

      if( i != Input -> Data.Nsim-1 ) {
	equate( &temp , Input -> Data.y[j] ) ;
	divide( &temp , Input -> Data.y[j+Input->Data.Ndata[i]] ) ;
      } else {
	equate( &temp , Input -> Data.y[j-Input->Data.Ndata[i]] ) ;
	divide( &temp , Input -> Data.y[j] ) ;
      }
      res_log( &temp ) ;

      if( Input -> Fit.Fitdef != CORNELL_V2 ) {
	equate( &Input -> Data.y[j] , temp ) ;
      }
      
      fprintf( effmass , "%e %e %e\n" ,
	       Input->Data.x[j].avg , temp.avg , temp.err ) ;
      
    }
    fprintf( effmass , "\n" ) ;
    shift += Input -> Data.Ndata[i] ;
  }
  fclose( effmass ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;  
  
  if( Input -> Fit.Fitdef == POLY ) {
    root( &Fit[1] ) ;
    fprintf( stdout , "rsigma %e %e\n" , Fit[1].avg , Fit[1].err ) ;
  } else if( Input -> Fit.Fitdef == CORNELL ) {
    root( &Fit[0] ) ;
    fprintf( stdout , "rsigma %e %e\n" , Fit[0].avg , Fit[0].err ) ;
  } else if( Input -> Fit.Fitdef == CORNELL_V2 ) {
    root( &Fit[1] ) ;
    fprintf( stdout , "rsigma %e %e\n" , Fit[1].avg , Fit[1].err ) ;
  }

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
