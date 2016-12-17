#include "gens.h"

#include "stats.h"
#include "bootfit.h"
#include "correlation.h"

int main( void )
{
  gsl_rng_env_setup( ) ;
  gsl_rng *r = gsl_rng_alloc( gsl_rng_default ) ;

  const size_t n = 40 , NBOOTS = 500 ;

  struct resampled *x = malloc( n * sizeof( struct resampled ) ) ;
  struct resampled *y = malloc( n * sizeof( struct resampled ) ) ;

  size_t i , j ;
  for( i = 0 ; i < n ; i++ ) {

    x[i].resampled = malloc( NBOOTS * sizeof( double ) ) ;
    x[i].restype   = BOOTDATA ;
    x[i].NSAMPLES  = NBOOTS ;

    y[i].resampled = malloc( NBOOTS * sizeof( double ) ) ;
    y[i].restype   = BOOTDATA ;
    y[i].NSAMPLES  = NBOOTS ;

    const double k = gsl_rng_uniform( r ) * n ;

    for( j = 0 ; j < NBOOTS ; j++ ) {
      x[i].resampled[j] = k + gsl_ran_gaussian( r , 0.01 ) ;
      const double yy = 5 * exp ( -0.1 * x[i].resampled[j] ) ;
      y[i].resampled[j] = yy * ( 1 + gsl_ran_gaussian( r , 0.05 ) ) ; //+ 1.0 ;
    }

    compute_err( &y[i] ) ;
    compute_err( &x[i] ) ;
    printf( "%g || %g %g \n" , x[i].avg , y[i].avg , y[i].err ) ;
  }

  const corrtype CORRFIT = UNCORRELATED ;
  const double **W = (const double**)correlations_inv( y , CORRFIT , n ) ;

  //write_corrmatrix( W , n ) ;

  struct resampled *fit = perform_bootfit( x , y , W , 
					   n , 64 , EXP , CORRFIT ) ;

  make_xmgrace_graph( "test.agr" , "x" , "y" ) ;

  plot_data( x , y , n ) ;

  plot_data( x , y , n ) ;

  plot_fitfunction( fit , EXP , x , n , 64 , CORRFIT ) ;

  close_xmgrace_graph( ) ;

  for( i = 0 ; i < n ; i++ ) {
    free( x[i].resampled ) ; free( y[i].resampled ) ; free( (void*)W[i] ) ;
  }
  free( x ) ; free( y ) ; free( (void*)W ) ;

  gsl_rng_free( r ) ;

  return 0 ;
}
