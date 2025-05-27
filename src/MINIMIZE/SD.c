/**
   line search steepest-descent does a fast job at doing badly

   You can improve it for an exponential/cosh by dividing the gradient associated with the mass by a large (1E6) constant. Why this is is that a small change in the mass makes a big change in the overall fit
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h"
#include "summation.h"

#include <assert.h>

#define VERBOSE

// steepest-descent iterations this sucks don't use it
int
sd_iter( void *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  //assert( !"Routine only exists for educational purposes" ) ;
  
  // point to the fit descriptor
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  
  // counters and max iterations SDMAX
  const size_t SDMAX = 80000 ;
  size_t iters = 0 ;

  // allocate the temporary fitfunction for computing new steps 
  // down descent direction in the line search
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , Fit -> f.N ) ;
  copy_ffunction( &f2 , Fit->f ) ;
  
  // allocate the gradient
  double grad[ Fit -> Nlogic ] ;
  
  // evaluate the function, its first and second derivatives
  f2.Prior = Fit -> f.Prior = Fit -> Prior ;
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
  
  double chisq_diff = 10 , alpha = 1 ;
  while( sqrt(chisq_diff) > TOL && iters < SDMAX ) {
    // compute the derivative of the \chi^2 function
    get_gradient( grad , W , Fit ) ;    
    // line search along it
    alpha = line_search( &f2 , Fit -> f , grad , grad , *Fit , data , W ) ;
    // update fparams
    for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
      Fit -> f.fparams[i] += alpha*grad[i] ;
    }
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
    Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
    const double chi = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
    chisq_diff = fabs( chi - Fit -> f.chisq ) ;
    Fit -> f.chisq = chi ;   
    iters++ ;
  }

#ifdef VERBOSE
  // tell us how many iterations we hit
  if( iters == SDMAX ) {
    printf( "\n[SD] stopped by max iterations %zu \n" , iters ) ;
  } else {
    printf( "\n[SD] FINISHED in %zu iterations \n" , iters ) ;
  }
  printf( "[SD] chisq :: %e -> DIFF %e \n\n" , Fit -> f.chisq , chisq_diff ) ;
  for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
    printf( "PARAMS :: %e \n" , Fit -> f.fparams[i] ) ;
  }
#endif

  // free the gradient and fitfunction
  free_ffunction( &f2 , f2.NPARAMS ) ;

  return iters ;
}
