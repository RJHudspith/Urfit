/**
   line search steepest-descent
   does a fast job at doing badly
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h"
#include "summation.h"

//#define VERBOSE

// steepest-descent iterations this sucks don't use it
int
sd_iter( void *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // point to the fit descriptor
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  
  // counters and max iterations SDMAX
  const size_t SDMAX = 100000 ;
  size_t iters = 0 , i ;

  // allocate the temporary fitfunction for computing new steps 
  // down descent direction in the line search
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , Fit -> f.N ) ;

  // allocate the gradient
  double *grad = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  
  // evaluate the function, its first and second derivatives
  f2.Prior = Fit -> f.Prior = Fit -> Prior ;
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
  
  double chisq_diff = 10 , alpha = 1 ;;
  while( sqrt(chisq_diff) > TOL && iters < SDMAX ) {
    // compute the derivative of the \chi^2 function
    get_gradient( grad , W , Fit ) ;
    // line search along it
    alpha = line_search( &f2 , Fit -> f , grad , grad ,
			 *Fit , data , W , 100*alpha ) ;
    // update fparams
    for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
      Fit -> f.fparams[i] += alpha*grad[i] ;
    }
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
    Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
    const double chi = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
    chisq_diff = fabs( chi - Fit -> f.chisq ) ;
    Fit -> f.chisq = chi ;
    #ifdef VERBOSE
    printf( "[SD] ITER %zu :: chidiff %e \n" , iters , chisq_diff ) ;
    #endif 
    iters++ ;
  }

  // tell us how many iterations we hit
  if( iters == SDMAX ) {
    printf( "\n[SD] stopped by max iterations %zu \n" , iters ) ;
  } else {
    printf( "\n[SD] FINISHED in %zu iterations \n" , iters ) ;
  }

  printf( "[SD] chisq :: %e -> DIFF %e \n\n" , Fit -> f.chisq , chisq_diff ) ;
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    printf( "PARAMS :: %e \n" , Fit -> f.fparams[i] ) ;
  }

  // free the gradient and fitfunction
  free( grad ) ;
  free_ffunction( &f2 , f2.NPARAMS ) ;

  return iters ;
}
