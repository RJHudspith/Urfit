/**
   @file line_search.c
   @brief perform a backtracking line search
 */
#include "gens.h"

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
  const double fac = 0.5 ;
  double abest = alpha * 1000 ;
  size_t iters = 0 , ITERMAX = 50 ;

  double *y = NULL , *t = NULL ;
  size_t Nsum = f1.N ;
  switch( f1.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = f1.N * f1.N ; break ;
  }
  y = malloc( Nsum * sizeof( double ) ) ;
  
  while( iters < ITERMAX ) {
    //compute sd step 
    const double chi = test_step( f2 , descent[jidx] , f1 , fdesc , 
				  data , W , jidx , abest ) ;

    // add the armijo condition
    // == alpha * beta * g * p , where p is the descent direction and g is 
    // gradient is g
    double armijo = 0.0 , curve1 = 0.0 , curve2 = 1.0 ;

    const double c1 = 0.0001 , c2 = 0.1 ;
    armijo = c1 * abest * grad[jidx] * descent[jidx] ;
    #ifdef VERBOSE
    printf( "[LINE SEARCH] %e %e %e\n" , chi , f1.chisq , armijo ) ;
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
    printf( "[LINE SEARCH] LS diff :: %e \n" ,
	    fabs( chi - ( f1.chisq + armijo ) ) ) ;
    #endif

    // compute test chi  
    if( ( chi <= ( f1.chisq + armijo ) ) && 
	( curve1 >= curve2 ) ) {
      break ;
    }

    iters++ ;
    abest *= fac ;
  }

  // free the temporary storage
  if( y != NULL ) {
    free( y ) ;
  }
  
#ifdef VERBOSE
  if( iters == ITERMAX ) {
    printf( "[LINE SEARCH] backtracking failed\n" ) ;
  }
  printf( "[LINE SEARCH] abest :: %e \n" , abest ) ;
#endif
  return abest ;
}
