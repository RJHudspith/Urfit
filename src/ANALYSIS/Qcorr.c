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

//#define SUBZERO

//#define ZRESCALE

// the the topological corr < Q >^2 / V
int
fit_Qcorr( struct input_params *Input )
{   
  size_t i , j , shift = 0 ;
  const double fac = 1.0/ ( 10*10*10*10. ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift = j ;
  }
  
  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}

// partially sum the correlator \sum_{x=0}^{r^2} < q(x+y) q(x) >
int
fit_Qsusc( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // subtract the zero by hand?
    #ifdef SUBZERO
    equate_constant( &Input -> Data.y[shift] , 0.0 ,
		     Input -> Data.y[shift].NSAMPLES ,
		     Input -> Data.y[shift].restype ) ;
    #endif
    
    // compute the partial sum up to whatever cutoff
    for( j = shift+1 ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      
      add( &Input -> Data.y[j] , Input -> Data.y[j-1] ) ;
    }
    
    shift = j ;
  }

  average_equivalent( Input ) ;

  // compute the derivative of the partial sum
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim; i++ ) {

    struct resampled der = init_dist( NULL ,
				      Input -> Data.x[shift].NSAMPLES ,
				      Input -> Data.x[shift].restype ) ;

    for( j = shift+1 ; j < shift + Input -> Data.Ndata[i] - 1 ; j++ ) {

	equate( &der , Input -> Data.y[j+1] ) ;
	subtract( &der , Input -> Data.y[j-1] ) ;
	mult_constant( &der , 1.0/( Input -> Data.x[j+1].avg - Input -> Data.x[j-1].avg ) ) ;
	fprintf( stdout , "[DER] %e %e %e \n" , Input -> Data.x[j].avg , der.avg , der.err ) ;

	if( ( der.avg - der.err ) < 0.0 ) {
	  fprintf( stdout , "CANDIDATE %e %e %e \n" , Input -> Data.x[j].avg ,
		   Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
	  break ;
	}
    }

    free( der.resampled ) ;
    
    shift = j ;
  }

  // multiply
  shift = 0 ;
  const double fac = 1.0 / ( 10*10*10*10. ) ;
  for( i = 0 ; i < Input -> Data.Nsim; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift = j ;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
	
  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}

// fit the slab method evaluation of Qsusc
int
fit_Qslab( struct input_params *Input )
{
  size_t i = 0 , j , shift = 0 ;

  // V factor
  struct resampled *t0 = read_flat_single( "b17.074NC.flat" ) ;

  printf( "bootstrapping t0\n" ) ;
  
  bootstrap_single( &t0[0] , Input -> Data.Nboots ) ;
  //raise( &t0[0] , 2 ) ;
  const double t0sq = 1/(t0[0].avg) ;
  raise( &t0[0] , 2 ) ;

  fprintf( stdout , "t0^4 read %e %e\n" , t0[0].avg , t0[0].err ) ;

#ifdef ZRESCALE
  struct resampled *Z = read_flat_single( "HYP120b6.0945.flat" ) ;
  printf( "bootstrapping Z\n" ) ;
  bootstrap_single( &Z[0] , Input -> Data.Nboots ) ;
  raise( &Z[0] , 2 ) ;
  printf( "Z^2 %e %e\n" , Z[0].avg , Z[0].err ) ;
  mult( &t0[0] , Z[0] ) ;
#endif
  
  const double V = ( Input -> Traj[0].Dimensions[0] *
		     Input -> Traj[0].Dimensions[1] *
		     Input -> Traj[0].Dimensions[2] *
		     Input -> Traj[0].Dimensions[3] ) ;

  divide_constant( &t0[0] , V ) ;
  mult_constant( &t0[0] , (double)Input -> Traj[0].Dimensions[3] ) ;

  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    // SU2 guy
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      add_constant( &Input -> Data.x[j] , 1.0 ) ;
      mult( &Input -> Data.y[j] , t0[0] ) ;
    }
    shift=j;
  }

  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    struct resampled res = init_dist( NULL ,
				      Input -> Data.y[shift].NSAMPLES ,
				      Input -> Data.y[shift].restype ) ;
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      if( j == (shift + Input -> Data.Ndata[i] - 1) ) continue ;
      equate( &res , Input -> Data.y[j+1] ) ;
      subtract( &res , Input -> Data.y[j] ) ;
      printf( "%e %e %e \n" , Input -> Data.x[j].avg ,
	      res.avg , res.err ) ;
    }
    shift=j;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  if( Input -> Fit.Fitdef == QSLAB ) {
    //mult_constant( &fit[2] , Input -> Traj[0].Dimensions[3] / fac ) ;
    //printf( "T0 * meta :: %e %e \n" , fit[2].avg , fit[2].err ) ;
  }

  // compute the usual Q^2
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    divide_constant( &Input -> Data.y[j] , Input -> Traj[0].Dimensions[3] ) ;
  }
  printf( "<Q2> %e %e\n" ,
	  Input->Data.y[ Input->Data.Ndata[0]-1 ].avg ,
	  Input->Data.y[ Input->Data.Ndata[0]-1 ].err ) ;

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

// clean these up
#ifdef SUBZERO
  #undef SUBZERO
#endif
