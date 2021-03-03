/**
   @file Qcorr.c
   @brief topological correlator measurement
 */
#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "momenta.h"
#include "read_flat.h"
#include "resampled_ops.h"
#include "bootstrap.h"

// fit the slab method evaluation of Qsusc
int
fit_Qslab_fixed( struct input_params *Input )
{
  size_t i = 0 , j , shift = 0 ;

  // V factor
  //struct resampled *t0 = read_flat_single( "../b24.786NC.flat" ) ;
  //struct resampled *t0 = read_flat_single( "../b25.175NC.flat" ) ;
  //struct resampled *t0 = read_flat_single( "b45.74NC.flat" ) ;
  struct resampled *t0 = read_flat_single( "b57.369NC.flat" ) ;
  //struct resampled *t0 = read_flat_single( "b70.682NC.flat" ) ;
  //struct resampled *t0 = read_flat_single( "b71.25NC.flat" ) ;
  //struct resampled *t0 = read_flat_single( "b101.788NC.flat" ) ;

  printf( "bootstrapping t0\n" ) ;
  
  bootstrap_single( &t0[0] , Input -> Data.Nboots ) ;
  const double t0sq = 1/t0[0].avg ;
  raise( &t0[0] , 2 ) ;

  const double Q = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // V factor
    const double V = 
      (const double)(Input -> Traj[i].Dimensions[0]*
		     Input -> Traj[i].Dimensions[1]*
		     Input -> Traj[i].Dimensions[2]*
		     Input -> Traj[i].Dimensions[3]) ;
    
    // Factor of 4 is because I didn't remove the ND in the computation
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      add_constant( &Input -> Data.x[j] , 1+Q*1. ) ; // is this right? Looks it
      // rescale x to the interval (0,1]
      divide_constant( &Input -> Data.x[j] ,
		       (double)Input -> Traj[i].Dimensions[3] ) ;
      mult_constant( &Input -> Data.y[j] , 4 ) ;
      const double x = Input -> Data.x[j].avg ;
      subtract_constant( &Input -> Data.y[j] , Q*x*x ) ;
      // normalise by V so the result of fitparam[0] is chi - susceptibility
      mult_constant( &Input -> Data.y[j] , 1./V ) ;
      mult( &Input -> Data.y[j] , t0[0] ) ;
    }
    shift=j;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  // write out a flat file?
  char str[256] ;
  sprintf( str , "Qslab_L%zu.flat" , Input -> Traj[0].Dimensions[0] ) ;
  FILE *file = fopen( str , "w" ) ;
  fprintf( file , "%u\n" , fit[0].restype ) ;
  fprintf( file , "1\n" ) ;
  fprintf( file , "%zu\n" , fit[0].NSAMPLES ) ;
  for( j = 0 ; j < Input -> Data.y[0].NSAMPLES ; j++ ) {
    fprintf( file , "%1.15e %1.15e\n" , t0sq , fit[0].resampled[j] ) ;
  }
  fprintf( file , "AVG %1.15e %1.15e\n" , t0sq , fit[0].avg ) ;
  fclose( file ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}
