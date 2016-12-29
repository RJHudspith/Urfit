/**
   @file Mainfile.c
   @brief Mainly a mainfile
 */
#include "gens.h"

#include "bootfit.h"
#include "correlation.h"
#include "make_xmgrace.h"
#include "plot_fitfunc.h"
#include "pmap.h"
#include "stats.h"

int main( void )
{
  gsl_rng_env_setup( ) ;
  gsl_rng *r = gsl_rng_alloc( gsl_rng_default ) ;

  const size_t Nsim      = 4 ;
  const size_t Nboots    = 500 ;
  const size_t n[ 4 ]    = { 10 , 10 , 10 , 10 } ;
  const double amp[ 4 ]  = { 200 , 200 , 200 , 200 } ;
  const double mass[ 4 ] = { 1.0 , 1.0 , 1.0 , 1.0 } ;
  bool S[ 2 ] = { true , true } ;
  
  size_t h , Ndata = 0  ;
  for( h = 0 ; h < Nsim ; h++ ) {
    Ndata += n[h] ;
  }

  struct resampled *x = malloc( Ndata * sizeof( struct resampled ) ) ;
  struct resampled *y = malloc( Ndata * sizeof( struct resampled ) ) ;

  size_t i , j , shift = 0 ;
  for( h = 0 ; h < Nsim ; h++ ) {
    for( i = shift ; i < shift + n[h] ; i++ ) {

      const size_t idx = i ;
      x[idx].resampled = malloc( Nboots * sizeof( double ) ) ;
      x[idx].restype   = BOOTDATA ;
      x[idx].NSAMPLES  = Nboots ;

      y[idx].resampled = malloc( Nboots * sizeof( double ) ) ;
      y[idx].restype   = BOOTDATA ;
      y[idx].NSAMPLES  = Nboots ;

      const double k = gsl_rng_uniform( r ) * n[h] ;

      for( j = 0 ; j < Nboots ; j++ ) {
	x[idx].resampled[j] = k + gsl_ran_gaussian( r , 0.05 ) ;
	const double yy = amp[h] * exp ( -mass[h] * x[idx].resampled[j] ) ;
	y[idx].resampled[j] = yy * ( 1 + gsl_ran_gaussian( r , 0.05 ) ) ;
      }

      compute_err( &y[idx] ) ;
      compute_err( &x[idx] ) ;
      printf( "%g || %g %g \n" , x[idx].avg , y[idx].avg , y[idx].err ) ;
    }
    shift += n[h] ;
  }

  const corrtype CORRFIT = UNCORRELATED ;
  const double **W = (const double**)correlations_inv( y , CORRFIT , Ndata ) ;
  struct resampled *fit = NULL ;

  write_corrmatrix( (const double**)W , Ndata ) ;

  size_t Npars ;
  fit = perform_bootfit( &Npars , x , y , W , n , Nsim ,
			 S , 64 , EXP , CORRFIT ) ;

  // make the graph
  make_xmgrace_graph( "test.agr" , "x" , "y" ) ;

  shift = 0 ;
  for( h = 0 ; h < Nsim ; h++ ) {
    plot_data( x + shift , y + shift , n[h] ) ;
    shift += n[h] ;
  }
  
  plot_fitfunction( fit , EXP , x , n ,
		    64 , CORRFIT , Nsim , S ) ;

  close_xmgrace_graph( ) ;

  // free all the data
  for( i = 0 ; i < Ndata ; i++ ) {
    free( x[i].resampled ) ; free( y[i].resampled ) ; free( (void*)W[i] ) ;
  }
  free( x ) ; free( y ) ; free( (void*)W ) ;

  // free the fit
  for( i = 0 ; i < Npars ; i++ ) {
    free( fit[i].resampled ) ;
  }
  free( fit ) ;
  
  gsl_rng_free( r ) ;

  return 0 ;
}
