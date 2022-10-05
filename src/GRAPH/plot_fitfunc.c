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
#include "resampled_ops.h" // compute err
#include "stats.h"

#include "fvol2.h"
#include "fvol3.h"

struct resampled
extrap_fitfunc( const struct resampled *f ,
		const struct data_info Data ,
		const struct fit_info Fit ,
		const double xpos ,
		const size_t shift )
{
  // set up the fit again
  struct fit_descriptor fdesc = init_fit( Data , Fit ) ;

  struct resampled data = init_dist( NULL , f[0].NSAMPLES , f[0].restype ) ;
  
  struct x_desc xdesc = { xpos , Data.LT[shift] , Fit.N , Fit.M } ;
  //struct x_desc xdesc = { xpos , 0.0 , Fit.N , Fit.M } ;
  double fparams[ fdesc.Nparam ] ;
  
  size_t j , p ;
  for( j = 0 ; j < f[0].NSAMPLES ; j++ ) {
    // evaluate the fitfunc
    for( p = 0 ; p < fdesc.Nparam ; p++ ) {
      fparams[ p ] = f[ Fit.map[shift].p[p] ].resampled[j] ;
      //fparams[p] = f[ p ].resampled[j] ;
    }
    data.resampled[j] = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;
    //data.resampled[j] = fdesc.func( xdesc , fparams , shift ) ;
  }
  // and the average
  for( p = 0 ; p < fdesc.Nparam ; p++ ) {
    fparams[ p ] = f[ Fit.map[shift].p[p] ].avg ;
    //fparams[ p ] = f[ p ].avg ;
  }
  data.avg = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;
  //data.avg = fdesc.func( xdesc , fparams , shift ) ;
  compute_err( &data ) ;
	
  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  return data ;
}

double
bisect_fitfunc( const struct resampled *f ,
		const struct data_info Data ,
		const struct fit_info Fit ,
		const double xlo ,
		const double target ,
		const double xhi ,
		const size_t shift )
{
  struct resampled tavg = init_dist( NULL , f[0].NSAMPLES , f[0].restype ) ;
  double x1 = xlo , x2 = xhi , tol = 1 , xmid = 0.0 ;
  while( tol > 1E-5 ) {
    xmid = 0.5*(x1+x2) ;
    tavg = extrap_fitfunc( f , Data , Fit , xmid , 0 ) ;
    if( tavg.avg > target ) {
      x2 = xmid ;
    } else {
      x1 = xmid ;
    }
    tol = fabs(x2-x1) ;
    //fprintf( stdout , "(%f,%f,%f) ----> %f\n", target, tavg.avg, tol, xmid ) ;
  }
  free( tavg.resampled ) ;
  return xmid ;
}

static struct resampled
extrap_fitfunc_HACK( const struct resampled *f ,
		     const struct data_info Data ,
		     const struct fit_info Fit ,
		     const double xpos ,
		     const size_t shift )
{
  // set up the fit again
  struct fit_descriptor fdesc = init_fit( Data , Fit ) ;

  struct resampled data = init_dist( NULL , f[0].NSAMPLES , f[0].restype ) ;

  struct x_desc xdesc = { xpos , 0.0 , Fit.N , Fit.M } ;
  double fparams[ fdesc.Nparam ] ;
  
  size_t j , p ;
  for( j = 0 ; j < f[0].NSAMPLES ; j++ ) {
    // evaluate the fitfunc
    for( p = 0 ; p < fdesc.Nparam ; p++ ) {
      fparams[p] = f[ p ].resampled[j] ;
    }
    data.resampled[j] = fdesc.func( xdesc , fparams , shift ) ;
  }
  // and the average
  for( p = 0 ; p < fdesc.Nparam ; p++ ) {
    fparams[ p ] = f[ p ].avg ;
  }
  data.avg = fdesc.func( xdesc , fparams , shift ) ;
  compute_err( &data ) ;
	
  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  return data ;
}


int
plot_fitfunction_HACK( const struct resampled *f ,
		       const struct data_info Data ,
		       const struct fit_info Fit )
{    
  // loop x with this granularity
  const size_t granularity = 501 ;
  double *X    = malloc( granularity * sizeof( double ) ) ;
  double *YAVG = malloc( granularity * sizeof( double ) ) ;
  double *YMIN = malloc( granularity * sizeof( double ) ) ;
  double *YMAX = malloc( granularity * sizeof( double ) ) ;

  size_t h , i ;

  const size_t MAX = 13 ; //fvol2_NMAX( ) ;
  //const size_t MAX = 15 ; //fvol3_NMAX( ) ;

  // loop over the simultaneous parameters
  for( size_t shift = 0 ; shift < MAX ; shift++ ) {
    //  size_t shift = 0 ;
    for( h = 0 ; h < 1 ; h++ ) {
      
      //const double xmin = 0.018225 ;
      //const double xmax = 0.18 ;

      const double xmin = 0.07822077018599871;
      const double xmax = 0.8;
      
      //const double xmin = 0.225981 ;
      //const double xmax = 0.0 ;
      //const double xmin = 0.245025 ;
      //const double xmax = 0.16 ;
      
      const double x_step = ( xmax - xmin ) / ( granularity - 1 ) ;
      for( i = 0 ; i < granularity ; i++ ) {
	
	X[ i ] = xmin + x_step*i ;
	
	struct resampled data = extrap_fitfunc_HACK( f , Data , Fit , X[i] , shift ) ;
	
	YMAX[i] = data.err_hi ;
	YAVG[i] = data.avg ;
	YMIN[i] = data.err_lo ;
	
	// free the fit distribution
	free( data.resampled ) ;
      }

      // draw lines between the evaluated fitfunctions
      draw_line( X , YMAX , granularity ) ;
      draw_line( X , YAVG , granularity ) ;
      draw_line( X , YMIN , granularity ) ;

      if( shift == MAX-1 ) {
	fprintf( stdout , "Extrap result %e %e\n" , YAVG[0] , 0.5*(YMAX[0]-YMIN[0]) ) ;
      }
      
    }
  }

  // free the x, y , ymin and ymax
  free( X ) ; free( YAVG ) ; free( YMIN ) ; free( YMAX ) ;

  return SUCCESS ;
}

int
plot_feffmass( const struct resampled *f ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  // loop x with this granularity
  const size_t granularity = 501 ;
  double *X    = malloc( granularity * sizeof( double ) ) ;
  double *YAVG = malloc( granularity * sizeof( double ) ) ;
  double *YMIN = malloc( granularity * sizeof( double ) ) ;
  double *YMAX = malloc( granularity * sizeof( double ) ) ;

  size_t h , i ;
  FILE *file = fopen( "feffmass.dat" , "w" ) ;

  // loop over the simultaneous parameters
  size_t shift = 0 ;
  for( h = 0 ; h < Data.Nsim ; h++ ) {

    // loop the x to find the max and min of x
    double xmin = Data.x[shift].err_lo ;
    double xmax = Data.x[shift].err_hi ;
    for( i = shift ; i < shift + Data.Ndata[h] ; i++ ) {
      if( Data.x[i].err_lo < xmin ) {
	xmin = Data.x[i].err_lo ;
      }
      if( Data.x[i].err_hi > xmax ) {
	xmax = Data.x[i].err_hi ;
      }
    }

    const double x_step = ( xmax - xmin ) / ( granularity - 1 ) ;
    for( i = 0 ; i < granularity ; i++ ) {

      const double dh = 1E-5 ;
      
      X[ i ] = xmin + x_step*i ;
      const double X2 = xmin + x_step*i + dh ;

      // do a derivative here
      
      struct resampled data1 = extrap_fitfunc( f , Data , Fit , X[i] , shift ) ;
      struct resampled data2 = extrap_fitfunc( f , Data , Fit , X2 , shift ) ;

      divide( &data1 , data2 ) ;
      res_log( &data1 ) ;
      divide_constant( &data1 , dh ) ;

      YMAX[i] = data1.err_hi ;
      YAVG[i] = data1.avg ;
      YMIN[i] = data1.err_lo ;

      // free the fit distribution
      free( data1.resampled ) ;
      free( data2.resampled ) ;
    }

    // draw lines between the evaluated fitfunctions
    for( size_t i = 0 ; i < granularity ; i++ ) {
      fprintf( file , "%e %e\n" , X[i] , YMAX[i] ) ;
    }
    fprintf( file , "\n" ) ;
    for( size_t i = 0 ; i < granularity ; i++ ) {
      fprintf( file , "%e %e\n" , X[i] , YAVG[i] ) ;
    }
    fprintf( file , "\n" ) ;
    for( size_t i = 0 ; i < granularity ; i++ ) {
      fprintf( file , "%e %e\n" , X[i] , YMIN[i] ) ;
    }
    fprintf( file , "\n" ) ;
    
    shift += Data.Ndata[h] ;
  }
  fclose( file ) ;

  // free the x, y , ymin and ymax
  free( X ) ; free( YAVG ) ; free( YMIN ) ; free( YMAX ) ;

  return SUCCESS ;
}

int
plot_fitfunction( const struct resampled *f ,
		  const struct data_info Data ,
		  const struct fit_info Fit )
{    
  // loop x with this granularity
  const size_t granularity = 501 ;
  double *X    = malloc( granularity * sizeof( double ) ) ;
  double *YAVG = malloc( granularity * sizeof( double ) ) ;
  double *YMIN = malloc( granularity * sizeof( double ) ) ;
  double *YMAX = malloc( granularity * sizeof( double ) ) ;

  size_t h , i ;

  // loop over the simultaneous parameters
  size_t shift = 0 ;
  for( h = 0 ; h < Data.Nsim ; h++ ) {

    // loop the x to find the max and min of x
    double xmin = Data.x[shift].err_lo ;
    double xmax = Data.x[shift].err_hi ;
    for( i = shift ; i < shift + Data.Ndata[h] ; i++ ) {
      if( Data.x[i].err_lo < xmin ) {
	xmin = Data.x[i].err_lo ;
      }
      if( Data.x[i].err_hi > xmax ) {
	xmax = Data.x[i].err_hi ;
      }
    }
    //xmin = 0.11233596051497033 ;
    //xmin = 0 ;

    const double x_step = ( xmax - xmin ) / ( granularity - 1 ) ;
    for( i = 0 ; i < granularity ; i++ ) {
      
      X[ i ] = xmin + x_step*i ;
	    
      struct resampled data = extrap_fitfunc( f , Data , Fit , X[i] , shift ) ;

      YMAX[i] = data.err_hi ;
      YAVG[i] = data.avg ;
      YMIN[i] = data.err_lo ;

      if( i == 0 ) {
	fprintf( stdout , "XMIN %e %e\n" , data.avg , data.err ) ;
      }

      // free the fit distribution
      free( data.resampled ) ;
    }

    // draw lines between the evaluated fitfunctions
    draw_line( X , YMAX , granularity ) ;
    draw_line( X , YAVG , granularity ) ;
    draw_line( X , YMIN , granularity ) ;

    shift += Data.Ndata[h] ;
  }

  // free the x, y , ymin and ymax
  free( X ) ; free( YAVG ) ; free( YMIN ) ; free( YMAX ) ;

  return SUCCESS ;
}
