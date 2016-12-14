/**
   line_search.c
   perform a backtracking line search
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"

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
  // perform a "backtracking line search"
  double t = 0.6180339887498547 ;
  double abest = alpha ;
  size_t iters = 0 , ITERMAX = 50 ;
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
    printf( "%e %e %e\n" , chi , f1.chisq , armijo ) ;
    #endif

    size_t i , k ;
    register double newgrad = 0.0 ;
    for( i = 0 ; i < f1.N ; i++ ) {
      switch( f1.CORRFIT ) {
      case UNWEIGHTED :
	newgrad += -f2 -> df[jidx][i] * f2 -> f[i] ;
	break ;
      case UNCORRELATED :
	newgrad += -f2 -> df[jidx][i] * W[i][i] * f2 -> f[i] ;
	break ;
      case CORRELATED :
	for( k = 0 ; k < f1.N ; k++ ) {
	  newgrad += -f2 -> df[jidx][i] * W[i][k] * f2 -> f[k] ;
	}
	break ;
      }
    }
    // add any prior stuff
    if( f1.prior[ jidx ] != UNINIT_FLAG ) {
      newgrad -= ( f2 -> fparams[ jidx ] - f2 -> prior[ jidx ] ) /
	( f2 -> err_prior[ jidx ] * f2 -> err_prior[ jidx ] ) ;
    }
    curve1 = descent[jidx] * newgrad ;
    curve2 = c2 * descent[ jidx ] * grad[ jidx ] ;

    #ifdef VERBOSE
    printf( "LS diff :: %e \n" , fabs( chi - ( f1.chisq + armijo ) ) ) ;
    #endif

    // compute test chi  
    if( ( chi <= ( f1.chisq + armijo ) ) && 
	( curve1 >= curve2 ) ) {
      break ;
    }

    iters++ ;
    abest *= t ;
  }
#ifdef VERBOSE
  if( iters == ITERMAX ) {
    printf( "[CG] backtracking failed\n" ) ;
  }
  printf( "abest :: %e \n" , abest ) ;
#endif
  return abest ;
}
