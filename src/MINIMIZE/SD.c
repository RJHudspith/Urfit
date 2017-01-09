/**
   line search steepest-descent
   does a fast job at doing badly
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h" 

// largest step size we use
#define BIG_GUESS (10)

// steepest-descent iterations
int
sd_iter( struct fit_descriptor *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // counters and max iterations SDMAX
  const size_t SDMAX = 10000 ;
  size_t iters = 0 , i ;

  // allocate the temporary fitfunction for computing new steps 
  // down descent direction in the line search
  struct ffunction f2 = allocate_ffunction( fdesc -> Nlogic , 
					    fdesc -> f.N ) ;

  // allocate the gradient
  double *grad = malloc( fdesc -> Nlogic * sizeof( double ) ) ;

  // set the guesses
  fdesc -> guesses( fdesc -> f.fparams , fdesc -> Nlogic ) ;

  // get priors
  fdesc -> set_priors( fdesc -> f.prior , fdesc -> f.err_prior ) ;

  // evaluate the function, its first and second derivatives
  fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ;
  fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;
  fdesc -> f.chisq = compute_chisq( fdesc -> f , W , fdesc -> f.CORRFIT ) ;

  double chisq_diff = 10 ;
  while( chisq_diff > TOL && iters < SDMAX ) {

    // compute the descent direction ( - the gradient for SD! )
    for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
      grad[i] = 0.0 ;
      size_t j , k ;
      for( j = 0 ; j < fdesc -> f.N ; j++ ) {
	switch( fdesc -> f.CORRFIT ) {
	case UNWEIGHTED :
	  grad[i] += -fdesc -> f.df[i][j] * fdesc -> f.f[j] ;
	  break ;
	case UNCORRELATED :
	  grad[i] += -fdesc -> f.df[i][j] * W[j][j] * fdesc -> f.f[j] ;
	  break ;
	case CORRELATED :
	  for( k = 0 ; k < fdesc -> f.N ; k++ ) {
	    grad[i] += -fdesc -> f.df[i][j] * W[j][k] * fdesc -> f.f[k] ;
	  }
	  break ;
	}
      }
      if( fdesc -> f.prior[i] != UNINIT_FLAG ) {
	grad[i] -= ( fdesc -> f.fparams[i] - fdesc -> f.prior[i] ) / 
	  ( fdesc -> f.err_prior[i] * fdesc -> f.err_prior[i] ) ;
      }
    }
    for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
      // line search best alpha
      const double ap = line_search( &f2 , fdesc -> f , grad , grad ,
				     *fdesc , data , W , i , BIG_GUESS ) ;
      fdesc -> f.fparams[i] += ap * grad[i] ;

      #ifdef VERBOSE
      printf( "[SD] line ap :: %e \n" , ap ) ;
      #endif
    }

    fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ;
    fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;

    const double chi = compute_chisq( fdesc -> f , W , fdesc -> f.CORRFIT ) ;
    chisq_diff = fabs( chi - fdesc -> f.chisq ) ;
    fdesc -> f.chisq = chi ;

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

  printf( "[SD] chisq :: %e -> DIFF %e \n\n" , fdesc -> f.chisq , chisq_diff ) ;
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    printf( "PARAMS :: %f \n" , fdesc -> f.fparams[i] ) ;
  }

  // free the gradient
  free( grad ) ;

  // free the fitfunction
  free_ffunction( &f2 , f2.NPARAMS ) ;

  return iters ;
}

#undef BIG_GUESS
