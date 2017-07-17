/**
   @brief quick check of the exceptional fits
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "resampled_ops.h"

#include "init.h"

//#define DOUBLE_POLE
#define SINGLE_POLE
//#define LINEAR_FIT

int
fit_exceptional( struct input_params *Input )
{
  #ifdef LINEAR_FIT
  if( Input -> Fit.Fitdef != POLY && Input -> Fit.N == 1 ) {
    fprintf( stderr , "Linear fit expects linear polynomial fit\n" ) ;
    return FAILURE ;
  }
  #endif
  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      // convert am to physical units
      if( i%2 == 1 ) {
	mult_constant( &Input -> Data.x[j] , 2.383 ) ;
      } else {
	mult_constant( &Input -> Data.x[j] , 1.785 ) ;
      }
      
      #ifdef DOUBLE_POLE
      mult_constant( &Input -> Data.y[j] ,
		     Input -> Data.x[j].avg * Input -> Data.x[j].avg ) ;
      #elif ( defined SINGLE_POLE ) || (defined  LINEAR_FIT)
      mult_constant( &Input -> Data.y[j] , Input -> Data.x[j].avg ) ;
      #endif
      
    }
    shift = j ;
  }
  //exit(1) ;

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  #ifdef LINEAR_FIT
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    // subtract the pole
    subtract( &Input -> Data.y[i] , fit[0] ) ;
    // divide data through by m again
    mult_constant( &Input -> Data.y[i] , 1.0 / ( Input -> Data.x[i].avg ) ) ;
  }
  #endif

  // write out a flat file
  FILE *outfile1 = fopen( "Zinv.dat" , "w" ) ;
  
  fprintf( outfile1 , "%d\n" , fit[0].restype ) ;
  fprintf( outfile1 , "%zu\n" , 1 ) ;
  
  // write out a flat file
  fprintf( outfile1 , "%zu\n" , fit[0].NSAMPLES ) ;
  for( i = 0 ; i < fit[0].NSAMPLES ; i++ ) {
    #if (defined SINGLE_POLE)
    if( Input -> Fit.Fitdef == POLES ) {
      fprintf( outfile1 , "%1.12e %1.12e\n" , 0.0 , fit[5].resampled[i] ) ;
    } else {
      fprintf( outfile1 , "%1.12e %1.12e\n" , 0.0 , fit[1].resampled[i] ) ;
    }
    #elif (defined DOUBLE_POLE)
    fprintf( outfile1 , "%1.12e %1.12e\n" , 0.0 , fit[2].resampled[i] ) ;
    #elif !(defined LINEAR_FIT)
    fprintf( outfile1 , "%1.12e %1.12e\n" , 0.0 , fit[2].resampled[i] ) ;
    #endif
  }

  if( Input -> Fit.Fitdef == POLES ) {
    FILE *outfile2 = fopen( "Zinv2.dat" , "w" ) ;
    fprintf( outfile2 , "%d\n" , fit[0].restype ) ;
    fprintf( outfile2 , "%zu\n" , 1 ) ;
    fprintf( outfile2 , "%zu\n" , fit[0].NSAMPLES ) ;
    for( i = 0 ; i < fit[0].NSAMPLES ; i++ ) {
      fprintf( outfile2 , "%1.12e %1.12e\n" , 0.0 , fit[6].resampled[i] ) ;
    }
    fclose( outfile2 ) ;
  }

  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  #ifdef LINEAR_FIT

  Input -> Traj[0].Fit_High = 0.0425 ;
  
  fit = fit_and_plot( *Input , &chisq ) ;
  for( i = 0 ; i < fit[0].NSAMPLES ; i++ ) {
    fprintf( outfile1 , "%1.12e %1.12e\n" , 0.0 , fit[0].resampled[i] ) ;
  }
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  #endif

  fclose( outfile1 ) ;
  
  return SUCCESS ;
}
