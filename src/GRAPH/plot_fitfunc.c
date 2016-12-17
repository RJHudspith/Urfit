/**
   @file plot_fitfunc.c
   @brief plot the fitfuntion
 */
#include "gens.h"

#include "fit_chooser.h"
#include "fitfunc.h"

#include "make_xmgrace.h" // drawing lines

void
plot_fitfunction( const struct resampled *f ,
		  const fittype fit ,
		  const struct resampled *x ,
		  const size_t Ndata ,
		  const size_t LT ,
		  const corrtype CORRFIT )
{
  struct fit_descriptor fdesc = init_fit( fit , Ndata , CORRFIT ) ;

  // loop the x to find the max and min of x
  double xmin = x[0].err_lo , xmax = x[0].err_hi ;
  size_t i , j ;
  for( i = 0 ; i < Ndata ; i++ ) {
    if( x[i].resampled[j] < xmin ) {
      xmin = x[i].err_lo ;
    }
    if( x[i].resampled[j] > xmax ) {
      xmax = x[i].err_hi ;
    }
  }

  // loop x
  const size_t granularity = 101 ;
  double *X    = malloc( granularity * sizeof( double ) ) ;
  double *YAVG = malloc( granularity * sizeof( double ) ) ;
  double *YMIN = malloc( granularity * sizeof( double ) ) ;
  double *YMAX = malloc( granularity * sizeof( double ) ) ;

  const double x_step = ( xmax - xmin ) / ( granularity - 1 ) ;
  for( i = 0 ; i < granularity ; i++ ) {
    const double extrap_x = xmin + x_step*i ;
    X[ i ] = extrap_x ;
    // evaluate the fitfunc
    struct x_desc xdesc = { X[i] , LT } ;
    double fparams[ fdesc.NPARAMS ] ;
    size_t p ;
    // compute the hi values
    for( p = 0 ; p < fdesc.NPARAMS ; p++ ) {
      fparams[ p ] = f[p].err_hi ;
    }
    YMAX[i] = fdesc.func( xdesc , fparams , fdesc.NPARAMS ) ;
    // compute the lo values
    for( p = 0 ; p < fdesc.NPARAMS ; p++ ) {
      fparams[ p ] = f[p].err_lo ;
    }
    YMIN[i] = fdesc.func( xdesc , fparams , fdesc.NPARAMS ) ;
    // compute the mid values
    for( p = 0 ; p < fdesc.NPARAMS ; p++ ) {
      fparams[ p ] = f[p].avg ;
    }
    YAVG[i] = fdesc.func( xdesc , fparams , fdesc.NPARAMS ) ;
  }

  // draw lines between the evaluated fitfunctions
  draw_line( X , YMAX , granularity ) ;
  draw_line( X , YAVG , granularity ) ;
  draw_line( X , YMIN , granularity ) ;

  // free the x, y , ymin and ymax
  free( X ) ; free( YAVG ) ; free( YMIN ) ; free( YMAX ) ;

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.NPARAMS ) ;

  return ;
}
