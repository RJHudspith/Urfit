/**
   @file plot_fitfunc.c
   @brief plot the fitfuntion
 */
#include "gens.h"

#include "fit_chooser.h"
#include "ffunction.h"
#include "fitfunc.h"
#include "make_xmgrace.h" // drawing lines
#include "pmap.h"

void
plot_fitfunction( const struct resampled *f ,
		  const fittype fit ,
		  const struct resampled *x ,
		  const size_t *Ndata ,
		  const size_t LT ,
		  const corrtype CORRFIT ,
		  const size_t Nsims ,
		  const bool *sims )
{
  size_t h , i , j , p , shift = 0 ;
  for( h = 0 ; h < Nsims ; h++ ) {
    shift += Ndata[h] ;
  }
    
  // set up the fit again
  struct fit_descriptor fdesc = init_fit( fit , shift , CORRFIT ,
					  Nsims , sims ) ;

  // set up the param map
  struct pmap *map = parammap( fdesc.Nparam , Nsims , Ndata , sims ) ;

  // loop x
  const size_t granularity = 101 ;
  double *X    = malloc( granularity * sizeof( double ) ) ;
  double *YAVG = malloc( granularity * sizeof( double ) ) ;
  double *YMIN = malloc( granularity * sizeof( double ) ) ;
  double *YMAX = malloc( granularity * sizeof( double ) ) ;

  // loop over the simultaneous parameters
  shift = 0 ;
  for( h = 0 ; h < Nsims ; h++ ) {

    // loop the x to find the max and min of x
    double xmin = x[shift].err_lo , xmax = x[shift].err_hi ;
    for( i = shift ; i < shift + Ndata[h] ; i++ ) {
      for( j = 0 ; j < x[i].NSAMPLES ; j++ ) {
	if( x[i].resampled[j] < xmin ) {
	  xmin = x[i].err_lo ;
	}
	if( x[i].resampled[j] > xmax ) {
	  xmax = x[i].err_hi ;
	}
      }
    }

    const double x_step = ( xmax - xmin ) / ( granularity - 1 ) ;
    for( i = 0 ; i < granularity ; i++ ) {
      const double extrap_x = xmin + x_step*i ;
      X[ i ] = extrap_x ;
      // evaluate the fitfunc
      struct x_desc xdesc = { X[i] , LT } ;
      double fparams[ fdesc.Nparam ] ;
      // compute the hi values
      for( p = 0 ; p < fdesc.Nparam ; p++ ) {
	fparams[ p ] = f[ map[shift].p[p] ].err_hi ;
      }
      YMAX[i] = fdesc.func( xdesc , fparams , fdesc.Nparam ) ;
      // compute the lo values
      for( p = 0 ; p < fdesc.Nparam ; p++ ) {
	fparams[ p ] = f[ map[shift].p[p] ].err_lo ;
      }
      YMIN[i] = fdesc.func( xdesc , fparams , fdesc.Nparam ) ;
      // compute the mid values
      for( p = 0 ; p < fdesc.Nparam ; p++ ) {
	fparams[ p ] = f[ map[shift].p[p] ].avg ;
      }
      YAVG[i] = fdesc.func( xdesc , fparams , fdesc.Nparam ) ;
    }

    // draw lines between the evaluated fitfunctions
    draw_line( X , YMAX , granularity ) ;
    draw_line( X , YAVG , granularity ) ;
    draw_line( X , YMIN , granularity ) ;

    shift += Ndata[h] ;
  }

  // free the x, y , ymin and ymax
  free( X ) ; free( YAVG ) ; free( YMIN ) ; free( YMAX ) ;

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  // free the parameter map
  free_pmap( map , shift ) ;

  return ;
}
