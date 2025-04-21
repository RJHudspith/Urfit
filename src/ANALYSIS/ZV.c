/**
   @file ZV.c
   @brief compute ZV from a constant fit to cl/ll vector correlators

   ll is actually <V>^2 / Z_V^2
   cl is acutally <V>^2 / Z_V

   so cl/ll should plateau to Z_V
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"

#include "write_flat.h"

//#define OBC
//#define OBCTLEN (96)

// takes the ratio of two correlators
int
ZV_old_analysis( struct input_params *Input )
{
  if( Input -> Data.Nsim != 2 ) {
    return FAILURE ;
  }
  
  // square the first data point
  size_t j ;
  const size_t N = Input -> Data.Ndata[0] ;

  for( j = 0 ; j < N ; j++ ) {
    divide( &Input -> Data.y[ j ] , Input -> Data.y[ j+N ] ) ;
  }

  Input -> Data.Nsim = 1 ;
  Input -> Data.Ntot = Input->Data.Ndata[0] ;
  Input -> Fit.N = Input -> Fit.M = 1 ;
  Input -> Fit.Nlogic = 2 ;

  // perform a fit to a constant
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}

// 3pt divided by the pion at t/2
int
ZV_analysis( struct input_params *Input )
{
  if( Input -> Data.Nsim != 2 ) {
    return FAILURE ;
  }
  
  // square the first data point
  size_t j ;
  const size_t T1 = Input -> Data.Ndata[0] ;
  #ifdef OBC
  const size_t T2 = OBCTLEN+1 ; //64+1 ; //48+1 ;
  #else
  const size_t T2 = Input -> Data.Ndata[1] ;
  #endif
  
  for( j = 0 ; j < T1 ; j++ ) {
    divide( &Input -> Data.y[j] , Input -> Data.y[T1+T2-1] ) ;
    raise( &Input -> Data.y[j] , -1 ) ;
    #ifdef OBC
    if( j > T2 ) {
      equate_constant( &Input->Data.y[j] , 0.0 ,
		       Input->Data.y[j].NSAMPLES ,
		       Input->Data.y[j].restype ) ;
    }
    #endif
  }
  
  // perform a fit to a constant
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  write_flat_dist( Fit , Input->Data.x , 1 , "ZV.flat" ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}
