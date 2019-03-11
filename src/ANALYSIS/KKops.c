/**
   @file KKops.c
   @brief code computes the ratio O_i / O_1
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "resampled_ops.h"
#include "init.h"
#include "effmass.h"

#define SYMMETRIZE

#define SYMMETRIZE2

int
fit_ratios( struct input_params *Input )
{  
  size_t i , j ;
  const size_t LT = Input -> Traj[0].Dimensions[3] ;
  const size_t SRC = Input -> Traj[0].Dimensions[0] ;

  for( j = 0 ; j < LT ; j++ ) {
    divide( &Input -> Data.y[j] , Input -> Data.y[j+LT] ) ;
    if( Input -> Data.Nsim == 4 ) {
      equate( &Input -> Data.y[j+LT] , Input -> Data.y[j+2*LT] ) ;
      divide( &Input -> Data.y[j+LT] , Input -> Data.y[j+3*LT] ) ;
    }
  }

  // symmetrize?
  #ifdef SYMMETRIZE
  for( j = 0 ; j < LT/2 ; j++ ) {
    add( &Input -> Data.y[j] , Input -> Data.y[j+LT/2] ) ;
    divide_constant( &Input -> Data.y[j] , 2.0 ) ;
    equate_constant( &Input -> Data.y[j+LT/2] , 0.0 ,
		      Input -> Data.y[j+LT/2].NSAMPLES ,
		      Input -> Data.y[j+LT/2].restype ) ;
    if( Input -> Data.Nsim == 4 ) {
      add( &Input -> Data.y[j+LT] , Input -> Data.y[j+LT/2+LT] ) ;
      divide_constant( &Input -> Data.y[j+LT] , 2.0 ) ;
      equate_constant( &Input -> Data.y[j+LT/2+LT] , 0.0 ,
		       Input -> Data.y[j+LT/2+LT].NSAMPLES ,
		       Input -> Data.y[j+LT/2+LT].restype ) ;
    }
  }
#ifdef SYMMETRIZE2
  for( j = 1 ; j < LT/2 ; j++ ) {
    if( j > SRC/2 ) {
      equate_constant( &Input -> Data.y[j] , 0.0 ,
		       Input -> Data.y[j].NSAMPLES ,
		       Input -> Data.y[j].restype ) ;
    } else {
      add( &Input -> Data.y[j] , Input -> Data.y[(SRC-j)] ) ;
      divide_constant( &Input -> Data.y[j] , 2.0 ) ;
    }
    
    if( Input -> Data.Nsim == 4 ) {
      if( j > SRC/2 ) {
	equate_constant( &Input -> Data.y[j+LT] , 0.0 ,
			 Input -> Data.y[j+LT].NSAMPLES ,
			 Input -> Data.y[j+LT].restype ) ;
      } else {
	add( &Input -> Data.y[j+LT] , Input -> Data.y[(SRC-j)+LT] ) ;
	divide_constant( &Input -> Data.y[j+LT] , 2.0 ) ;
      }
    }
  }
  #endif
#endif

  Input -> Data.Nsim = Input -> Data.Nsim/2 ;
  Input -> Data.Ntot = Input -> Data.Ntot/2 ;
  
  struct resampled *effmass = effective_mass( Input , LOGBWD_EFFMASS ) ;
  
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;
  
  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  // write out the fitparam to a flat file
  char str[ 256 ] ;
  sprintf( str , "Q_%zu.flat" , Input -> Traj[0].Gs ) ;
  FILE *file = fopen( str , "w" ) ;
  fprintf( file , "%d\n" , Fit[0].restype ) ;
  fprintf( file , "%d\n" , 1 ) ;
  fprintf( file , "%zu\n" , Fit[0].NSAMPLES ) ;
  for( j = 0 ; j < Fit[0].NSAMPLES ; j++ ) {
    fprintf( file , "%1.15e %1.15e\n" , 0.0 , Fit[0].resampled[j] ) ;
  }
  fclose(file) ;
  

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  Input -> Data.Nsim = Input -> Data.Nsim*2 ;
  Input -> Data.Ntot = Input -> Data.Ntot*2 ;
  
  return SUCCESS ;
}

