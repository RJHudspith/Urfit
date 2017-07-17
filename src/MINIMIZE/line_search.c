/**
   @file line_search.c
   @brief perform a backtracking line search with Wolfe conditions
 */
#include "gens.h"

#include <limits.h>

#include "chisq.h"
#include "ffunction.h"
#include "summation.h"

//#define VERBOSE

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
  // really and probably will do at some point
  double fac = 0.1 ;

  // do a rough search 1^12 -> 1E-12 for abest step 100
  double min = 123456789 ;
  double atrial = 1E22 , abest = 1 ;
  while( atrial > 1E-15 ) {
    atrial *= 0.01 ;
    double trial = test_step( f2 , descent[jidx] , f1 , fdesc , 
			      data , W , jidx , atrial ) ;
    //printf( "%e trial %e min %e \n" , atrial , trial , min ) ;
    if( trial < min ) {
      abest = atrial ;
      min = trial ;
    }
  }

  size_t iters = 0 , ITERMAX = 25 ;

  double *y = NULL , *t = NULL ;
  size_t Nsum = f1.N ;
  switch( f1.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = f1.N * f1.N ; break ;
  }
  y = malloc( Nsum * sizeof( double ) ) ;

  double chia , chib , chic ;
  
  double amid = abest ;
  double aup = amid / fac , adn = amid * fac ;
  
  while( iters < ITERMAX ) {

    // bracket a minimum
    chia = test_step( f2 , descent[jidx] , f1 , fdesc , 
		      data , W , jidx , aup ) ;

    chib = test_step( f2 , descent[jidx] , f1 , fdesc , 
		      data , W , jidx , amid ) ;

    chic = test_step( f2 , descent[jidx] , f1 , fdesc , 
		      data , W , jidx , adn ) ;

    if( isnan( chia ) || isnan( chib ) || isnan( chic ) ) {
      aup *= fac ; amid *= fac ; adn *= fac ; continue ;
    } else if( isinf( chia ) || isinf( chib ) || isinf( chic ) ) {
      aup *= fac ; amid *= fac ; adn *= fac ; continue ;
    }

    if( chia < chib && chib < chic ) {
      adn = amid ; amid = aup ; aup = aup / ( fac ) ; 
    } else if( chic < chib && chib < chia ) {
      aup = amid ; amid = adn ; adn = adn * fac ;
    } else {
      fac = sqrt( fac ) ;
      aup = amid/fac ; adn = amid*fac ;
    }

    abest = amid ;

    // add the armijo condition
    // == alpha * beta * g * p , where p is the descent direction and g is 
    // gradient is g
    double armijo = 0.0 , curve1 = 0.0 , curve2 = 1.0 ;

    const double c1 = 1.0E-4 , c2 = 0.1 ;
    armijo = c1 * abest * grad[jidx] * descent[jidx] ;
    #ifdef VERBOSE
    printf( "[LINE SEARCH] %e %e %e\n" , chib , f1.chisq , armijo ) ;
    #endif

    size_t i , k ;
    t = y ;
    switch( f1.CORRFIT ) {
    case UNWEIGHTED :
      for( i = 0 ; i < f1.N ; i++ ) {
        *t = -f2 -> df[jidx][i] * f2 -> f[i] ; t++ ;
      }
      break ;
    case UNCORRELATED :
      for( i = 0 ; i < f1.N ; i++ ) {
	*t = -f2 -> df[jidx][i] * W[0][i] * f2 -> f[i] ; t++ ;
      }
      break ;
    case CORRELATED :
      for( i = 0 ; i < f1.N ; i++ ) {
	for( k = 0 ; k < f1.N ; k++ ) {
	  *t = -f2 -> df[jidx][i] * W[i][k] * f2 -> f[k] ; t++ ;
	}
      }
      break ;
    }
    register double newgrad = kahan_summation( y , Nsum ) ;
    // add any prior stuff
    if( f1.Prior[ jidx ].Initialised == true ) {
      newgrad -= ( f2 -> fparams[ jidx ] - f2 -> Prior[ jidx ].Val ) /
	( f2 -> Prior[ jidx ].Err * f2 -> Prior[ jidx ].Err) ;
    }
    curve1 = descent[jidx] * newgrad ;
    curve2 = c2 * descent[ jidx ] * grad[ jidx ] ;

    #ifdef VERBOSE
    printf( "[LINE SEARCH] LS diff %zu :: %e %e \n" ,
	    iters , chib , ( f1.chisq + armijo ) ) ;
    printf( "[LINE SEARCH] curves  %e %e \n" ,
	    curve1 , curve2 ) ;
    #endif

    // compute test chi  
    if( ( chib <= ( f1.chisq + armijo ) )
	//&& ( curve1 >= curve2 ) curve condition is tricky to get to work
	) {
      break ;
    }

    
    // some little checks to avoid wasting time
    if( abest < 1E-15 ) iters = ITERMAX ;

    iters++ ;
  }

  // free the temporary storage
  if( y != NULL ) {
    free( y ) ;
  }
  
  #ifdef VERBOSE
  if( iters == ITERMAX ) {
    fprintf( stderr , "[LINE SEARCH] backtracking failed %e \n" , abest ) ;
  }
  printf( "[LINE SEARCH] abest :: %e \n" , abest ) ;
  #endif

  return abest ;
}
