/**
   @file fit_and_plot.c
   @brief fit and plot the data after filtering it into the correct fit range
 */
#include "gens.h"

#include "bootfit.h"
#include "correlation.h"
#include "decays.h"
#include "init.h"         // free fit
#include "make_xmgrace.h" // drawing graph
#include "plot_fitfunc.h" 
#include "pmap.h"
#include "resampled_ops.h"
#include "stats.h"

#include "GLS.h"

// filter the data to lie within the fit range
static bool *
filter( size_t *N ,
	const struct data_info Data ,
	const struct fit_info Fit ,
	const struct traj *Traj )
{
  bool *in_fitrange = malloc( Data.Ntot * sizeof( bool ) ) ;
  size_t idx = 0 , i , j ;

  *N = 0 ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    for( j = 0 ; j < Data.Ndata[i] ; j++ ) {
      in_fitrange[ idx ] = false ;
      if( Data.x[ idx ].err_lo >= Traj[i].Fit_Low &&
	  Data.x[ idx ].err_hi <= Traj[i].Fit_High ) {
	in_fitrange[ idx ] = true ;
	*N = *N + 1 ;
      }
      idx ++ ;
    }
  }
  return in_fitrange ;
}

// performs a fit and plots the result
struct resampled *
fit_and_plot( struct input_params Input ,
	      double *Chi )
{
  // create a copy of the data struct
  struct data_info Data ;
  struct resampled *fitparams = NULL ;
  bool *in_fitrange = NULL ;

  Data.Nsim = Input.Data.Nsim ;
  in_fitrange = filter( &Data.Ntot , Input.Data ,
			Input.Fit , Input.Traj ) ;

  size_t i , j , idx = 0 ;

  Data.Cov.W = NULL ;
  Data.Restype = Input.Data.Restype ;
  Data.Cov.Eigenvalue_Tol = Input.Data.Cov.Eigenvalue_Tol ;
  Data.Cov.Column_Balanced = Input.Data.Cov.Column_Balanced ;
  Data.Cov.Divided_Covariance = Input.Data.Cov.Divided_Covariance ;

  Data.Ndata = NULL ;
  Data.x = NULL ;
  Data.y = NULL ;
  Data.LT = NULL ;

  // sanity check
  if( Data.Ntot == 0 ) {
    fprintf( stderr , "[FIT AND PLOT] Ntot is zero\n" ) ;
    goto memfree ;
  }
  
  Data.Ndata = malloc( Data.Nsim * sizeof( size_t ) ) ;
  Data.Nboots = Input.Data.Nboots ;

  Data.x = malloc( Data.Ntot * sizeof( struct resampled ) ) ;
  Data.y = malloc( Data.Ntot * sizeof( struct resampled ) ) ;
  
  size_t k = 0 ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    size_t Ndata = 0 ;
    for( j = 0 ; j < Input.Data.Ndata[i] ; j++ ) {
      if( in_fitrange[ k ] == true ) {
	Data.x[idx].resampled = malloc( Input.Data.x[k].NSAMPLES * sizeof( double ) ) ;
	Data.y[idx].resampled = malloc( Input.Data.x[k].NSAMPLES * sizeof( double ) ) ;
	equate( &Data.x[idx] , Input.Data.x[k] ) ;
	equate( &Data.y[idx] , Input.Data.y[k] ) ;

	printf( "%f %f %f\n" , Data.x[idx].avg , Data.y[idx].avg , Data.y[idx].err ) ;
	
	Ndata++ ; idx++ ;
      }
      k ++ ;
    }
    Data.Ndata[i] = Ndata ;
  }
  
  // set Lt
  if( init_LT( &Data , Input.Traj ) == FAILURE ) {
    fprintf( stderr , "[FIT AND PLOT] LT initialisation failure\n" ) ;
    goto memfree ;
  }

  // if we have a param map allocated we remove it and redo it
  // with the new ranges included?
  if( Input.Fit.map != NULL ) {
    free_pmap( Input.Fit.map , Input.Data.Ntot ) ;
  }
  Input.Fit.map = parammap( Data , Input.Fit ) ;

  if( inverse_correlation( &Data , Input.Fit ) == FAILURE ) {
    goto memfree ;
  }

#ifdef VERBOSE
  write_corrmatrix( (const double**)Data.Cov.W ,
		    Data.Ntot , Input.Fit.Corrfit ) ;
#endif

  // If we don't specify a fit
  if( Input.Fit.Fitdef != NOFIT ) {
    if( ( fitparams = perform_bootfit( Data , Input.Fit , Chi ) ) == NULL ) {
      goto memfree ;
    }
  }
  
  // make the graph
  make_xmgrace_graph( Input.Graph.Name ,
		      Input.Graph.Xaxis ,
		      Input.Graph.Yaxis ) ;

  size_t shift = 0 ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    plot_data( Input.Data.x + shift ,
	       Input.Data.y + shift ,
	       Input.Data.Ndata[i] ) ;
    shift += Input.Data.Ndata[i] ;
  }

  if( fitparams != NULL ) {
    plot_fitfunction( fitparams , Data , Input.Fit ) ;
  }
  
  close_xmgrace_graph( ) ;

 memfree :

  // free the list saying what is in the fitrange
  if( in_fitrange != NULL ) {
    free( in_fitrange ) ;
  }
  
  // free the data?
  free_Data( &Data , Input.Fit ) ;
  
  // free the parameter map
  if( Input.Fit.map != NULL ) {
    free_pmap( Input.Fit.map , Data.Ntot ) ;
  }
  
  return fitparams ;
}
