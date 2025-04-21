/**
   @file line_search.c
   @brief perform a backtracking line search with Wolfe conditions
 */
#include "gens.h"

#include <limits.h>
#include <assert.h>

#include "chisq.h"
#include "ffunction.h"
#include "summation.h"

// perform an step in the descent direction
static double
test_step( struct ffunction *f2 , 
	   const double grad ,
	   const struct ffunction f1 ,
	   const struct fit_descriptor fdesc ,
	   const void *data ,
	   const double **W ,
	   const size_t jidx ,
	   const double alpha )
{
  // copy f1 to f2
  copy_ffunction( f2 , f1 ) ;
  f2 -> fparams[ jidx ] += alpha * grad ;
  fdesc.F( f2 -> f , data , f2 -> fparams ) ;
  fdesc.dF( f2 -> df , data , f2 -> fparams ) ;
  return compute_chisq( *f2 , W , f2 -> CORRFIT ) ;
}

// backtracking line search with Wolfe conditions
double
line_search( struct ffunction *f2 ,
	     const struct ffunction f1 ,
	     const double *grad ,
	     const double *descent , 
	     const struct fit_descriptor fdesc ,
	     const void *data ,
	     const double **W ,
	     const size_t jidx ,
	     const double alpha )
{
  // perform a "backtracking line search" should use brent's method
  // really and probably will do at some point, uses a golden ratio
  // search to find approximately the minimum as long as this minimum
  // satisfies the wolfe conditions
  const double fac = 0.1 ;

  // do a rough search 1^12 -> 1E-12 for abest step 100
  double min = 1234567891000 ;
  double atrial = 1E8 , abest = 1E-15 ;
  size_t iters = 0 ;
  while( atrial > (1E-50)*alpha ) {
    atrial *= fac ;
    double trial = test_step( f2 , descent[jidx] , f1 , fdesc , 
			      data , W , jidx , atrial ) ;
    iters++ ;
    if( isnan( trial ) ) continue ;
    if( isinf( trial ) ) continue ;
    if( trial < min ) {
      abest = atrial ;
      min = trial ;
    }
  }
  double amid = abest ;
  double aup = amid / fac , adn = amid * fac ;
  // perform golden ratio line search
  const double gr = ( sqrt( 5. ) + 1. ) / 2. ;
  double c = aup - ( aup - adn ) / gr ;
  double d = adn + ( aup - adn ) / gr ;
  while( fabs( (c - d)/adn ) > 5E-3 ) {
    const double fc = test_step( f2 , descent[jidx] , f1 , fdesc , data , W , jidx , c ) ;
    const double fd = test_step( f2 , descent[jidx] , f1 , fdesc , data , W , jidx , d ) ;
    if( fc < fd ) {
      aup = d ;
    } else {
      adn = c ;
    }
    c = aup - ( aup - adn ) / gr ;
    d = adn + ( aup - adn ) / gr ;
    abest = ( aup + adn ) / 2. ;
  }
  #ifdef VERBOSE
  if( iters == 1000 ) {
    fprintf( stderr , "[LINE SEARCH] backtracking failed %e \n" , abest ) ;
  }
  printf( "[LINE SEARCH] abest :: %e \n" , abest ) ;
  #endif
  return abest ;
}
