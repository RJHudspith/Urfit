/**
   @file adler.c
   @brief adler function analysis
 */
#include "gens.h"

#include "resampled_ops.h"
#include "fit_and_plot.h"

//available a^-1
// 2.454315559701493 -> Z_V = 0.9553
// 3.607440054844607 -> Z_V = 0.9636
// 4.494919612756264 -> Z_V = 0.9699

int
adler_analysis( struct input_params *Input )
{
  const double ainv[3] = { 2.454315559701493 , 3.607440054844607 , 4.494919612756264 } ;
  const double Z_V[3] = { 0.9553 , 0.9636 , 0.9699 } ;
  size_t i , j , shift = 0 ;
  //for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( i = 0 ; i < 1 ; i++ ) {
      for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      // sqrt the x
      root( &( Input -> Data.x[j] ) ) ;
      mult_constant( &( Input -> Data.x[j] ) , ainv[0] ) ;
      /*
      // mult by Z_V^2
      mult_constant( &( Input -> Data.y[j] ) , Z_V[0]*Z_V[0] ) ;
      // multiply by 4\pi^2
      mult_constant( &( Input -> Data.y[j] ) , (4.*M_PI*M_PI) ) ;
      // subtract 1
      subtract_constant( &( Input -> Data.y[j] ) , 1.0 ) ;
      */
      printf( "Subtracting %zu\n" , j ) ;
      subtract( &( Input -> Data.y[j] ) ,
		Input -> Data.y[j+Input -> Data.Ndata[i] ] ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
