#include "gens.h"

#include "ffunction.h"

#include "LM.h"
#include "SD.h"
#include "CG.h"
#include "GA.h"

#include "exp.h"

int main( void )
{
  const size_t n = 40 ;

  gsl_rng_env_setup( ) ;
  gsl_rng *r = gsl_rng_alloc( gsl_rng_default ) ;

  double y[ n ] ;
  struct data d = { n, y } ;

  double **W = malloc( n * sizeof( double* ) ) ;

  size_t i ;
  for( i = 0 ; i < n ; i++ ) {
    W[i] = malloc( n * sizeof( double ) ) ;

    double t = i;
    double yi = 1.0 + 5 * exp (-0.1 * t);
    double si = 0.1 * yi;
    double dy = gsl_ran_gaussian(r, si);

    size_t j ;
    for( j = 0 ; j < n ; j++ ) {
      W[i][j] = 1.0 / ( ( 1 + fabs( j - i ) ) * si * si);
    }
    y[i] = yi + dy;
    //printf ("data: %zu %g %g\n", i, y[i], si);
    printf( "{%zu,%f},\n" , i , y[i] ) ;
  }

  struct fit_descriptor fdesc ;
  fdesc.F = exp_f ;
  fdesc.dF = exp_df ;
  fdesc.d2F = exp_d2f ;
  fdesc.guesses = exp_guesses ; 
  fdesc.set_priors = exp_priors ;
  fdesc.NPARAMS = 3 ;

  fdesc.f = allocate_ffunction( fdesc.NPARAMS , n ) ;
  fdesc.f.CORRFIT = CORRELATED ;

  int (*f) ( struct fit_descriptor *fdesc ,
	     const void *data ,
	     const double **W ,
	     const double TOL ) ;
  f = ga_iter ;

  f( &fdesc , &d , (const double**)W , 1E-10 ) ;

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.NPARAMS ) ;

  for( i = 0 ; i < n ; i++ ) {
    free( W[i] ) ;
  }
  free( W ) ;

  gsl_rng_free( r ) ;

  return 0 ;
}
