/**
   @file hvp_pade.c
   @brief fit a pade to the hvp
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "fit_chooser.h"
#include "ffunction.h"
#include "init.h"
#include "make_xmgrace.h"
#include "resampled_ops.h"
#include "stats.h"

void
pade_derivative( struct resampled *adler ,
		 const double q2 , 
		 const struct resampled *fit ,
		 const struct input_params *Input ,
		 const size_t shift )
{
  // do the individual bootstraps
  size_t j ;
  
  for( j = 0 ; j < fit[0].NSAMPLES ; j++ ) {
    double denom = 0.0 , num = 0.0 ;
    size_t n , m ;
    // left hand side
    for( n = 0 ; n < Input -> Fit.N ; n++ ) {
      num += fit[ n ].resampled[j] * pow( q2 , n ) ;
    }
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      denom += fit[ Input -> Fit.N + m ].resampled[j] * pow( q2 , (m+1) ) ;
    }
    denom += 1.0 ;
    // left hand side
    double left_sum = 0.0 , right_sum = 0.0 ;
    for( n = 1 ; n < Input -> Fit.N ; n++ ) {
      left_sum += n * fit[ n ].resampled[j] * pow( q2 , n-1 ) ;
    }
    left_sum = left_sum / denom ;
    // right hand side
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      right_sum += (m+1) * fit[ Input -> Fit.N + m ].resampled[j] * pow( q2 , m ) ;
    }
    right_sum = right_sum * num / ( denom * denom ) ;

    //printf( "%f -> %f %f \n" , num/denom , left_sum , right_sum ) ;
    adler -> resampled[j] = left_sum - right_sum ;
  }

  // do the average
  {
    double denom = 0.0 , num = 0.0 ;
    size_t n,m ;
    // left hand side
    for( n = 0 ; n < Input -> Fit.N ; n++ ) {
      num += fit[ n ].avg * pow( q2 , n ) ;
    }
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      denom += fit[ Input -> Fit.N + m ].avg * pow( q2 , (m+1) ) ;
    }
    denom += 1.0 ;
    // left hand side
    double left_sum = 0.0 , right_sum = 0.0 ;
    for( n = 1 ; n < Input -> Fit.N ; n++ ) {
      left_sum += n * fit[ n ].avg * pow( q2 , n-1 ) ;
    }
    left_sum = left_sum / denom ;
    // right hand side
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      right_sum += (m+1) * fit[ Input -> Fit.N + m ].avg * pow( q2 , m ) ;
    }
    right_sum = right_sum * num / ( denom * denom ) ;

    //printf( "%f -> %f %f \n" , num/denom , left_sum , right_sum ) ;
    adler -> avg = left_sum - right_sum ;
  }
  
  compute_err( adler ) ;
  //mult_constant( adler , 12*M_PI*M_PI*5./9. ) ;
  
  return ;
}

// numerical derivative
void
pade_derivative2( struct resampled *adler ,
		  const double q2 , 
		  const struct resampled *fit ,
		  const struct input_params *Input ,
		  const size_t shift )
{
  // set up the fit again
  struct fit_descriptor fdesc = init_fit( Input -> Data , Input -> Fit ) ;
  double fparams[ fdesc.Nparam ] ;
  size_t j , p ;
  const double h = 1E-7 ;

  // x+h and x-h
  struct x_desc xph = { q2 + h , 0 , Input -> Fit.N , Input -> Fit.M } ;
  struct x_desc xmh = { q2 - h , 0 , Input -> Fit.N , Input -> Fit.M } ;

  for( j = 0 ; j < fit[0].NSAMPLES ; j++ ) {
    // evaluate the fitfunc
    for( p = 0 ; p < fdesc.Nparam ; p++ ) {
      fparams[ p ] = fit[p].resampled[j] ;
    }
    adler -> resampled[j] =
      ( fdesc.func( xph , fparams , 0 ) -
	fdesc.func( xmh , fparams , 0 ) ) / ( 2 * h ) ;
  }

  // and the average
  for( p = 0 ; p < fdesc.Nparam ; p++ ) {
    fparams[ p ] = fit[p].avg ;
  }
  adler -> avg = ( fdesc.func( xph , fparams , 0 ) -
		   fdesc.func( xmh , fparams , 0 ) ) / ( 2 * h ) ;

  compute_err( adler ) ;
  //mult_constant( adler , q2*12*M_PI*M_PI*5./9. ) ;

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  return ;
}

int
fit_hvp( struct input_params *Input )
{
  size_t i = 0 ;
  double chisq ;

  // take the log of both sides
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    //raise( &Input -> Data.x[i] , 0.5 ) ;
  }
  
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  //exit(1) ;
  
  // compute the derivative q^2 d/dq^2
  struct resampled adler ;
  adler = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;

  make_xmgrace_graph( "der.agr" , "Q\\S2\\N [GeV\\S2\\N]" , "D(Q\\S2\\N)" ) ;
  
  const double inc = 0.01 ;
  const double hi = 1.0 ;
  const double lo = 0.0 ;
  const size_t Ndata = (size_t)( ( hi - lo ) / inc + 0.5 ) ;
  double *Q = malloc( Ndata * sizeof( double ) ) ;
  double *Ave = malloc( Ndata * sizeof( double ) ) ;
  double *Hi = malloc( Ndata * sizeof( double ) ) ;
  double *Lo = malloc( Ndata * sizeof( double ) ) ;

  size_t idx = 0 ;
  double q2 ;
  for( q2 = lo ; q2 < hi ; q2 += inc ) {

    pade_derivative2( &adler , q2 , fit , Input , i*Input -> Fit.Nparam ) ;

    Q[idx] = q2 ;
    Ave[idx] = adler.avg ; Hi[idx] = adler.err_hi ; Lo[idx] = adler.err_lo ;
    
    printf( "ADLER :: %f %f %f \n" , q2 , adler.avg , adler.err ) ;
    idx++ ;
  }

  draw_line( Q , Hi , Ndata ) ;
  draw_line( Q , Ave , Ndata ) ;
  draw_line( Q , Lo , Ndata ) ;

  free( Q ) ; free( Ave ) ; free( Hi ) ; free( Lo ) ;

  free( adler.resampled ) ;
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
    
  return SUCCESS ;
}
