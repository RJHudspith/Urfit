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

// armijo condition
int
armijo( struct ffunction *f2 ,
	const struct ffunction f1 ,
	const double *grad ,
	const double *descent , 
	const struct fit_descriptor fdesc ,
	const void *data ,
	const double **W ,
	const size_t jidx ,
	const double abest )
{
  double *y = NULL , *t = NULL ;
  size_t Nsum = f1.N ;
  switch( f1.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = f1.N * f1.N ; break ;
  }
  y = malloc( Nsum * sizeof( double ) ) ;
  
  // add the armijo condition
  // == alpha * beta * g * p , where p is the descent direction and g is 
  // gradient is g
  double armijo = 0.0 , curve1 = 0.0 , curve2 = 1.0 ;

  const double c1 = 1.0E-4 , c2 = 0.1 ;
  armijo = c1 * abest * grad[jidx] * descent[jidx] ;

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

  const double chib = test_step( f2 , descent[jidx] , f1 , fdesc , 
				 data , W , jidx , abest ) ;
#ifdef VERBOSE
  printf( "[LINE SEARCH] LS diff :: %e %e \n" ,
	  chib , ( f1.chisq + armijo ) ) ;
  printf( "[LINE SEARCH] curves  %e %e \n" ,
	  curve1 , curve2 ) ;
#endif
  
  // free the temporary storage
  if( y != NULL ) {
    free( y ) ;
  }

  return ( chib <= ( f1.chisq + armijo ) ) && ( curve1 >= curve2 ) ;
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
  double atrial = 1E8*alpha , abest = 1 ;
  size_t iters = 0 ;
  while( atrial > (1E-8)*alpha ) {
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

  const double gr = ( sqrt( 5. ) + 1. ) / 2. ;
  
  // perform golden ratio line search
  double c = aup - ( aup - adn ) / gr ;
  double d = adn + ( aup - adn ) / gr ;
  
  while( fabs( (c - d) ) > 1E-4 ) {
    const double fc = test_step( f2 , descent[jidx] , f1 , fdesc , 
				 data , W , jidx , c ) ;
    const double fd = test_step( f2 , descent[jidx] , f1 , fdesc , 
				 data , W , jidx , d ) ;
    if( fc < fd ) {
      aup = d ;
    } else {
      adn = c ;
    }
    
    c = aup - ( aup - adn ) / gr ;
    d = adn + ( aup - adn ) / gr ;

    if( armijo( f2 , f1 , grad , descent , fdesc , data , W ,
		jidx , ( aup + adn )/2 ) ) {
      abest = ( aup + adn ) / 2. ;
    }
  }
  
  #ifdef VERBOSE
  if( iters == ITERMAX ) {
    fprintf( stderr , "[LINE SEARCH] backtracking failed %e \n" , abest ) ;
  }
  printf( "[LINE SEARCH] abest :: %e \n" , abest ) ;
  #endif

  return abest ;
}
