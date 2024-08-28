/**
   @file effmass.c
   @brief effective mass solutions
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "gens.h"
#include "fitfunc.h"
#include "make_xmgrace.h"
#include "stats.h"
#include "resampled_ops.h"

#define EFFMASSTOL 1E-30

// set it to zero
static void
zero_effmass( struct resampled *effmass ,
	      const struct resampled bootavg )
{
  equate_constant( effmass , 0.0 , 
		   bootavg.NSAMPLES , 
		   bootavg.restype ) ;
  return ;
}

// asinh breaks down if acosh( n < 1 )
static bool
acosh_range( const struct resampled effmass )
{
  size_t i ;
  for( i = 0 ; i < effmass.NSAMPLES ; i++ ) {
    if( effmass.resampled[ i ] < 1 ) {
      return true ;
    }
  }
  return false ;
}

// make sure this is in the right range for the atanh
static bool
atanh_range( const struct resampled effmass ) 
{
  size_t i ;
  for( i = 0 ; i < effmass.NSAMPLES ; i++ ) {
    if( effmass.resampled[ i ] < -1 ||
	effmass.resampled[ i ] > 1 ) {
      return true ;
    }
  }
  return false ;
}

// if the data doesn't work we set the whole thing to zero
static bool
log_range( const struct resampled log )
{
  size_t i ;
  for( i = 0 ; i < log.NSAMPLES ; i++ ) {
    if( log.resampled[ i ] < 0 ) {
      return true ;
    }
  }
  return false ;
}

// log effective mass
static void
log_effmass( struct resampled *effmass ,
	     struct resampled y1 , struct resampled y2 ,
	     struct resampled x1 , struct resampled x2 )
{
  if( y2.avg == 0.0 ) { zero_effmass( effmass , y2 ) ; return ; }
  if( y1.avg == 0.0 ) { zero_effmass( effmass , y1 ) ; return ; }
  
  equate( effmass , y1 ) ;
  divide( effmass , y2 ) ;
  
  // if it is still negative we set to zero
  if( log_range( *effmass ) == true ) {
    zero_effmass( effmass , y1 ) ;
  } else {
    res_log( effmass ) ;
  }

  struct resampled Delta = init_dist( &x1 , x1.NSAMPLES , x1.restype ) ;

  subtract( &Delta , x2 ) ;
  divide( effmass , Delta ) ;
  
  free( Delta.resampled ) ;
}

// atanh effmass
static void
atanh_effmass( struct resampled *effmass ,
	       struct resampled y1 ,
	       struct resampled y2 )
{
  // is E^{-m(t-1)} - E^{-m(t+1)}
  equate( effmass , y1 ) ;
  subtract( effmass , y2 ) ;
  
  // is E^{-m(t-1)} + E^{-m(t+1)}
  struct resampled temp = init_dist( &y1 , y1.NSAMPLES , y1.restype ) ;
  add( &temp , y2 ) ;
  
  // divide the two and take the tanh
  divide( effmass , temp ) ;

  // if it is still negative we set to zero
  if( atanh_range( *effmass ) == true ) {
    zero_effmass( effmass , y1 ) ;
  } else {
    res_atanh( effmass ) ;
  }

  //
  free( temp.resampled ) ;
  
  return ;
}

// acosh effmass
static void
acosh_effmass( struct resampled *effmass ,
	       struct resampled y1 ,
	       struct resampled y2 ,
	       struct resampled y3 )
{
  equate( effmass , y1 ) ;
  add( effmass , y3 ) ;
  divide( effmass , y2 ) ;
  mult_constant( effmass , 0.5 ) ;
    
  // if it is still negative we set to zero
  if( acosh_range( *effmass ) == true ) {
    zero_effmass( effmass , y2 ) ;
  } else {
    res_acosh( effmass ) ;
  }

  return ;
}

// asinh effective mass
static void
asinh_effmass( struct resampled *effmass ,
	       struct resampled y1 ,
	       struct resampled y2 ,
	       struct resampled y3 )
{
  if( y2.avg == 0.0 ) { zero_effmass( effmass , y2 ) ; return ; }
  
  // compute ( y[i+1] - y[i-1] / y[i] ) 
  equate( effmass , y1 ) ;
  subtract( effmass , y3 ) ;
  //add( effmass , y3 ) ;
  divide( effmass , y2 ) ;
  mult_constant( effmass , 0.5 ) ;
    
  // asinh is valid for all inputs apart from exact zeros
  res_asinh( effmass ) ;

  return ;
}

// newton's method to solve
// c(t)/c(t+1) - (e^{-meff*t}+e^{-meff*(T-t)})/(e^{-meff*(t+1)}+e^{-meff*(T-t-1)) == 0
static inline double
fcosh( const double meff , const double ct , const double ct1 , const double t , const double t1 , const double Lt )
{
  return ct/ct1 - (exp( -meff*t ) + exp( -meff*(Lt-t) ))/( exp( -meff*t1 ) + exp( -meff*(Lt-t1) ) ) ;
}

// newton's method to solve
// c(t)/c(t+1) - (e^{-meff*t}-e^{-meff*(T-t)})/(e^{-meff*(t+1)}-e^{-meff*(T-t-1)) == 0
static inline double
fsinh( const double meff , const double ct , const double ct1 , const double t , const double t1 , const double Lt )
{
  return ct/ct1 - (exp( -meff*t ) - exp( -meff*(Lt-t) ))/( exp( -meff*t1 ) - exp( -meff*(Lt-t1) ) ) ;
}

// newton iteration
static double
newton_meff( const double initialguess ,
	     const double ct ,
	     const double ct1 ,
	     const double t ,
	     const double t1 ,
	     const int Lt ,
	     double (*f)( const double meff , const double ct , const double ct1 , const double t , const double t1 , const double Lt ) )
{
  size_t iter = 0 ;
  double meff = initialguess , tol = 1 , dh = 1E-8 ;

  if( fabs(t-t1) < 1E-14 ) {
    fprintf( stderr , "Newton coshmeff will fail 1 %e ~= t1 %e\n" , t , t1 ) ;
    exit(1) ;
  }
  
  while( tol > 1E-8 && iter < 20 ) {
    const double num = f(meff,ct,ct1,t,t1,(double)Lt) ;
    const double den = ( f(meff+dh,ct,ct1,t,t1,(double)Lt)
			 - f(meff-dh,ct,ct1,t,t1,(double)Lt) )/(2.*dh);    
    meff -= num/(den) ;
    tol = fabs( num ) ;
    iter++ ;
  }
  if( isnan( meff ) || isinf( meff ) || iter == 20 ) {
    fprintf( stderr , "Newton coshmeff failed %zu %e %zu\n" , t , meff , iter ) ;
    //exit(1) ;
  }
  return meff ;
}

// iterative cosh
static void
iterative_mass( struct resampled *effmass ,
		struct resampled y1 , struct resampled y2 ,
		struct resampled x1 , struct resampled x2 ,
		const int LT ,
		double (*f)( const double meff , const double ct , const double ct1 , const double t , const double t1 , const double Lt ) )
{
  if( y2.avg == 0.0 ) { zero_effmass( effmass , y2 ) ; return ; }
  if( y1.avg == 0.0 ) { zero_effmass( effmass , y1 ) ; return ; }
  // do the avg first --> need a more educated guess
  effmass -> avg = newton_meff( 0.5 , y1.avg , y2.avg , x1.avg , x2.avg , LT , f ) ;
  for( size_t k = 0 ; k < effmass -> NSAMPLES ; k++ ) {
    effmass -> resampled[k] = newton_meff( effmass -> avg ,
					   y1.resampled[k] , y2.resampled[k] ,
					   x1.resampled[k] , x2.resampled[k] , LT ,
					   f ) ;
  }
  compute_err( effmass ) ;
}

// computes the effective mass and plots a graph of it
struct resampled *
effective_mass( struct input_params *Input ,
		const effmass_type type )
{
  make_xmgrace_graph( "effmass.agr" , "t/a" , "am\\seff" ) ;
  
  struct resampled *effmass = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;

  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      effmass[j].NSAMPLES = Input -> Data.y[j].NSAMPLES ;
      effmass[j].restype  = Input -> Data.y[j].restype ;
      effmass[j].resampled = malloc( Input -> Data.y[j].NSAMPLES * sizeof( double ) ) ;

      if( j == shift ) {
	log_effmass( &effmass[j] ,
		     Input -> Data.y[j] , Input -> Data.y[j+1] ,
		     Input -> Data.x[j+1] , Input -> Data.x[j] ) ;
      } else if( j == ( shift + Input -> Data.Ndata[i] - 1 ) ) {
	log_effmass( &effmass[j] ,
		     Input -> Data.y[j-1] , Input -> Data.y[j] ,
		     Input -> Data.x[j] , Input -> Data.x[j-1] ) ;
      } else {

	// switch the various effective masses
	switch( type ) {
	case LOG_EFFMASS :
	case LOGFWD_EFFMASS :
	  log_effmass( &effmass[j] ,
		       Input -> Data.y[j] , Input -> Data.y[j+1] ,
		       Input -> Data.x[j+1] , Input -> Data.x[j] ) ;
	  break ;
	case LOGBWD_EFFMASS :
	  log_effmass( &effmass[j] ,
		       Input -> Data.y[j-1] , Input -> Data.y[j] ,
		       Input -> Data.x[j-1] , Input -> Data.x[j] ) ;
	  mult_constant( &effmass[j] , -1 ) ;
	  break ;
	case ATANH_EFFMASS :
	  atanh_effmass( &effmass[j] ,
			 Input -> Data.y[j-1] , Input -> Data.y[j+1] ) ;
	  break ;
	case ACOSH_EFFMASS :
	  acosh_effmass( &effmass[j] ,
			 Input -> Data.y[j-1] ,
			 Input -> Data.y[j] ,
			 Input -> Data.y[j+1] ) ;
	  break ;
	case ACOSH_ITERATIVE_EFFMASS :
	  iterative_mass( &effmass[j] ,
			  Input -> Data.y[j] , Input -> Data.y[j+1] ,
			  Input -> Data.x[j] , Input -> Data.x[j+1] ,
			  Input -> Traj[i].Dimensions[3] , fcosh ) ;
	  break ;
	case ASINH_EFFMASS :
	  asinh_effmass( &effmass[j] ,
			 Input -> Data.y[j-1] ,
			 Input -> Data.y[j] ,
			 Input -> Data.y[j+1] ) ;
	  break ;
	case ASINH_ITERATIVE_EFFMASS :
	  iterative_mass( &effmass[j] ,
			  Input -> Data.y[j] , Input -> Data.y[j+1] ,
			  Input -> Data.x[j] , Input -> Data.x[j+1] ,
			  Input -> Traj[i].Dimensions[3] , fsinh ) ;
	  break ;
	default :
	  // do nothing
	  break ;
	}
      }
    }
    plot_data( Input -> Data.x + shift , effmass + shift ,
	       Input -> Data.Ndata[i] ) ;
    
    shift = j ;
  }

  close_xmgrace_graph() ;
  
  return effmass ;
}

#undef EFFMASSTOL
