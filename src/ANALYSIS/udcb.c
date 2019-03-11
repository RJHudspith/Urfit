/**
   @file udcb.c
   @brief fitting the heavy mass dependence of udcb
 */
#include "gens.h"

#include "init.h"
#include "fit_and_plot.h"
#include "resampled_ops.h"

#include "plot_fitfunc.h"
#include "pmap.h"

#define INVERSE

int
udcb_analysis( struct input_params *Input )
{
  size_t i ;
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    mult_constant( &(Input -> Data.y[i]) , -1.0 ) ;
  }
  
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  // this is really bad, please think about it!!
  printf( "Extrapolation\n" ) ;
  if( Input -> Fit.map == NULL ) {
    Input -> Fit.map = parammap( Input -> Data , Input -> Fit ) ;
  }

  double ex[ 4 ] = { 3.077 , 3.077 , 3.077 , 3.077  } ;
  size_t shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    struct resampled temp = extrap_fitfunc( Fit ,
					    Input -> Data ,
					    Input -> Fit ,
					    ex[i] ,
					    shift ) ;
    shift += Input -> Data.Ndata[i] ;
    printf( "Extrap_%zu (%e) %e %e\n" , i , ex[i] , temp.avg , temp.err ) ;
    free( temp.resampled ) ;
  }

  if( Input -> Fit.map != NULL ) {
    free_pmap( Input -> Fit.map , Input -> Data.Ntot ) ;
  }
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return 0 ;
}
