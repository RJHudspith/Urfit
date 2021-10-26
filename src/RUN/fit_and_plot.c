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
#include "write_flat.h"

#include "Nint.h"
#include "GLS.h"

#define NINTPAR

//
typedef enum { NO_ERROR , NTOT_FAIL ,
	       LT_FAIL , CORRELATION_FAIL } prune_errflag ;

// filter the data to lie within the fit range
static bool *
filter( size_t *N ,
	const struct data_info Data ,
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

// applies filters so that only the data within our fit range gets fit
static prune_errflag
prune_data_for_fit( struct input_params *Input ,
		    struct data_info *Data )
{
  bool *in_fitrange = NULL ;
  prune_errflag error = NO_ERROR ;
  
  Data -> Nsim = Input -> Data.Nsim ;
  in_fitrange = filter( &Data -> Ntot , Input -> Data , Input -> Traj ) ;

  size_t i , j , idx = 0 ;

  Data -> Cov.W = NULL ;
  Data -> Restype = Input -> Data.Restype ;
  Data -> Cov.Eigenvalue_Tol = Input -> Data.Cov.Eigenvalue_Tol ;
  Data -> Cov.Column_Balanced = Input -> Data.Cov.Column_Balanced ;
  Data -> Cov.Divided_Covariance = Input -> Data.Cov.Divided_Covariance ;

  Data -> Ndata = NULL ;
  Data -> x = NULL ;
  Data -> y = NULL ;
  Data -> LT = NULL ;

  // sanity check
  if( Data -> Ntot == 0 ) {
    fprintf( stderr , "[FIT AND PLOT] Ntot is zero\n" ) ;
    error = NTOT_FAIL ;
    goto memfree ;
  }
  
  Data -> Ndata  = malloc( Data -> Nsim * sizeof( size_t ) ) ;
  Data -> Nboots = Input -> Data.Nboots ;
  Data -> x = malloc( Data -> Ntot * sizeof( struct resampled ) ) ;
  Data -> y = malloc( Data -> Ntot * sizeof( struct resampled ) ) ;

  size_t k = 0 ;
  for( i = 0 ; i < Data -> Nsim ; i++ ) {
    size_t Ndata = 0 ;
    for( j = 0 ; j < Input -> Data.Ndata[i] ; j++ ) {
      if( in_fitrange[ k ] == true ) {
	Data -> x[idx].resampled =
	  malloc( Input -> Data.x[k].NSAMPLES * sizeof( double ) ) ;
	Data -> y[idx].resampled =
	  malloc( Input -> Data.x[k].NSAMPLES * sizeof( double ) ) ;

	equate( &Data -> x[idx] , Input -> Data.x[k] ) ;
	equate( &Data -> y[idx] , Input -> Data.y[k] ) ;
	//#ifdef verbose
	  fprintf( stdout , "%e %e %e\n" ,
		   Data -> x[idx].avg ,
		   Data -> y[idx].avg ,
		   Data -> y[idx].err ) ;
	//#endif
	Ndata++ ; idx++ ;
      }
      k ++ ;
    }
    Data -> Ndata[i] = Ndata ;
  }
  
  // set Lt
  if( init_LT( Data , Input -> Traj ) == FAILURE ) {
    fprintf( stderr , "[FIT AND PLOT] LT initialisation failure\n" ) ;
    error = LT_FAIL ;
    goto memfree ;
  }

  // if we have a param map allocated we remove it and redo it
  // with the new ranges included?
  if( Input -> Fit.map != NULL ) {
    free_pmap( Input -> Fit.map , Input -> Data.Ntot ) ;
  }
  Input -> Fit.map = parammap( *Data , Input -> Fit ) ;

  if( inverse_correlation( Data , Input -> Fit ) == FAILURE ) {
    error = CORRELATION_FAIL ;
    goto memfree ;
  }

#ifdef VERBOSE
  write_corrmatrix( (const double**)Data -> Cov.W ,
		    Data -> Ntot , Input -> Fit.Corrfit ) ;
#endif

 memfree :

  // free the list saying what is in the fitrange
  if( in_fitrange != NULL ) {
    free( in_fitrange ) ;
  }
  
  return error ;
}

// performs a fit and plots the result and returns the fit parameters
struct resampled *
fit_and_plot( struct input_params Input ,
	      double *Chi )
{
  // create a copy of the data struct
  prune_errflag error = NO_ERROR ;
  struct data_info Data ;
  error = prune_data_for_fit( &Input , &Data ) ;

  struct resampled *fitparams = NULL ;
  
  if( error != NO_ERROR ) {
    goto memfree ;
  }
  
  // If we don't specify a fit we can just continue with the plotting
  if( Input.Fit.Fitdef != NOFIT ) {
    if( ( fitparams = perform_bootfit( Data , Input.Fit , Chi ) ) == NULL ) {
      goto memfree ;
    }
  }

#if 0
  // su3
  double pos[ 4 ] = { 16.8762 , 17.1074 , 17.3386 , 17.5585 } ;
  const int LT[ 4 ] = { 5 , 6 , 7 , 8  } ;
  size_t x ;
  for( x = 0 ; x < 4 ; x++ ) {
    const double shft = pos[x]-17.432 ;
    
    struct resampled t0 = extrap_fitfunc( fitparams ,
					   Input.Data ,
					   Input.Fit ,
					   shft , 0 ) ;
    printf( "[T0] %f %f %f\n" , pos[x] , t0.avg , t0.err ) ;

    struct resampled a2 = init_dist( &t0 ,
				     t0.NSAMPLES ,
				     t0.restype ) ;
    raise( &a2 , -2 ) ;
    
    divide_constant( &t0 , LT[x] ) ;
    printf( "[TC] %f %f %f\n" , a2.avg , t0.avg , t0.err ) ;	     
    free( t0.resampled ) ;
  }
#endif
  
  // make the graph
  make_xmgrace_graph( Input.Graph.Name ,
		      Input.Graph.Xaxis ,
		      Input.Graph.Yaxis ) ;

  size_t shift = 0 , i ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    plot_data( Input.Data.x + shift ,
	       Input.Data.y + shift ,
	       Input.Data.Ndata[i] ) ;
    shift += Input.Data.Ndata[i] ;
  }

  if( fitparams != NULL ) {
    plot_fitfunction( fitparams , Data , Input.Fit ) ;
    //plot_fitfunction_HACK( fitparams , Data , Input.Fit ) ;

    plot_feffmass( fitparams , Data , Input.Fit ) ;
  }
  
  close_xmgrace_graph( ) ;

 memfree :
  
  // free the data?
  free_Data( &Data , Input.Fit ) ;
  
  // free the parameter map
  if( Input.Fit.map != NULL ) {
    free_pmap( Input.Fit.map , Data.Ntot ) ;
  }
  
  return fitparams ;
}

// performs a fit and plots the result and returns the fit parameters
struct resampled *
fit_and_plot_and_Nint( struct input_params Input ,
		       double *Chi )
{
  // create a copy of the data struct
  prune_errflag error = NO_ERROR ;
  struct data_info Data ;
  error = prune_data_for_fit( &Input , &Data ) ;

  struct resampled *fitparams = NULL ;
  
  if( error != NO_ERROR ) {
    goto memfree ;
  }
 
  // If we don't specify a fit we can just continue with the plotting
  if( Input.Fit.Fitdef != NOFIT ) {
    if( ( fitparams = perform_bootfit( Data , Input.Fit , Chi ) ) == NULL ) {
      goto memfree ;
    }
  }
  
  // make the graph
  make_xmgrace_graph( Input.Graph.Name ,
		      Input.Graph.Xaxis ,
		      Input.Graph.Yaxis ) ;

  size_t shift = 0 , i ;
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

#ifdef NINTPAR
  // numerically integrate fit parameters ?
  shift = 0 ;
  size_t idx = 0 , j ;

  // numerically integrate up to fit_low
  struct resampled *YL =
    malloc( (Input.Data.Ndata[0]+1)*sizeof( struct resampled ) ) ;
  // numerically integrate up to fit_low
  struct resampled *XL =
    malloc( (Input.Data.Ndata[0]+1)*sizeof( struct resampled ) ) ;

  YL[0] = init_dist( NULL ,
		     Input.Data.y[0].NSAMPLES ,
		     Input.Data.y[0].restype ) ;
  XL[0] = init_dist( NULL ,
		     Input.Data.y[0].NSAMPLES ,
		     Input.Data.y[0].restype ) ;
  
  size_t n ;
  for( n = 0 ; n < Input.Data.Ndata[0]+1 ; n++ ) {
    YL[n+1] = init_dist( &Input.Data.y[n] ,
			 Input.Data.y[n].NSAMPLES ,
			 Input.Data.y[n].restype ) ;
    XL[n+1] = init_dist( &Input.Data.x[n] ,
			 Input.Data.x[n].NSAMPLES ,
			 Input.Data.x[n].restype ) ;
    mult_constant( &YL[n+1] , pow( Input.Data.x[n].avg , 3 ) ) ;

    printf( "Lint test %e %e %e\n" , XL[n+1].avg , YL[n+1].avg , YL[n+1].err ) ;
    if( Input.Data.x[n].avg > Input.Traj[0].Fit_Low ) break ;
  }
  struct resampled LInt = Nint( XL , YL , n+2 , true ) ;
  for( size_t j = 0 ; j < n ; j++ ) {
    free( YL[j].resampled ) ;
    free( XL[j].resampled ) ;
  }
  free( XL ) ;
  free( YL ) ;

  fprintf( stdout, "CHECK LINT %e %e %e\n" ,
	   Input.Data.x[n].avg , LInt.avg , LInt.err ) ;

  double stp = 1 ;
  const size_t N = (size_t)((64-Data.x[0].avg)/stp)+1 ;
  
  struct resampled *Y = malloc( N*sizeof( struct resampled ) ) ;
  struct resampled *X = malloc( N*sizeof( struct resampled ) ) ;


  // and then do the fit for the rest
  for( size_t j = 0 ; j < n ; j++ ) {
    mult_constant( &Input.Data.y[j] ,
		   pow( Input.Data.x[j].avg , 3 ) ) ;
    printf( "Grand %e %e %e\n" , Input.Data.x[j].avg ,
	    Input.Data.y[j].avg , Input.Data.y[j].err ) ;
  }

  // and then do the fit for the rest
  idx=0 ;
  for( double x = Data.x[0].avg ; x < 64 ; x+=stp ) {
    Y[idx] = extrap_fitfunc( fitparams , Data ,
			     Input.Fit ,
			     x , 0 ) ;
    X[idx].resampled = malloc( Y[idx].NSAMPLES*sizeof(double) ) ;
    equate_constant( &X[idx] , x , Y[idx].NSAMPLES , Y[idx].restype ) ;

    //mult_constant( &Y[idx] , pow( x , 3 )/0.06426 ) ;
    mult_constant( &Y[idx] , pow( x , 3 ) ) ;
    printf( "Grand %e %e %e %e\n" , X[idx].avg , Y[idx].err_hi , Y[idx].avg , Y[idx].err_lo ) ;
    idx++ ;
  }
  
  for( n = 1 ; n < N ; n++ ) { 
    struct resampled Int = Nint( X , Y , n , true ) ;
    add( &Int , LInt ) ;
    fprintf( stdout , "Gral %e %e %e %e\n" , X[n-1].avg , Int.err_hi , Int.avg , Int.err_lo ) ;
    if( n == (N-1) ) {
      struct resampled mpi2 = init_dist( NULL ,
					 Int.NSAMPLES ,
					 Int.restype ) ;
      char str[256] ;
      sprintf( str , "Gral.flat" ) ;
      equate_constant( &mpi2 , 0.0 , Int.NSAMPLES , Int.restype ) ;
      write_flat_dist( &Int , &mpi2 , 1 , str ) ;
    }
    free( Int.resampled ) ;
  }
  for( i = 0 ; i < Input.Data.Nsim ; i++ ) {
    // loop the x to find the max and min of x as they are not sorted
    // we need to traverse the entire array
    double xmin = Data.x[shift].err_lo ;
    double xmax = Data.x[shift].err_hi ;
    for( j = shift ; j < shift + Data.Ndata[i] ; j++ ) {
      if( Data.x[j].err_lo < xmin ) {
	xmin = Data.x[j].err_lo ;
      }
      if( Data.x[j].err_hi > xmax ) {
	xmax = Data.x[j].err_hi ;
      }
    }
    struct resampled Int = Nint_fit( fitparams , Data , Input.Fit ,
				     xmax , xmin , 1E-12 , i ) ;
    free( Int.resampled ) ;
    shift += Data.Ndata[i] ;
  }

  free( LInt.resampled ) ;
  #endif

 memfree :
  
  // free the data?
  free_Data( &Data , Input.Fit ) ;
  
  // free the parameter map
  if( Input.Fit.map != NULL ) {
    free_pmap( Input.Fit.map , Data.Ntot ) ;
  }
  
  return fitparams ;
}
