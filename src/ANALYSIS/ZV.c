/**
   @file ZV.c
   @brief compute ZV from a constant fit to cl/ll vector correlators

   ll is actually <V>^2 / Z_V^2
   cl is acutally <V>^2 / Z_V

   so cl/ll should plateau to Z_V
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"

int
ZV_analysis( struct input_params *Input )
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

  //free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}
