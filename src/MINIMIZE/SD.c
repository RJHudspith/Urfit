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

// steepest-descent iterations
int
sd_iter( void *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // point to the fit descriptor
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  
  // counters and max iterations SDMAX
  const size_t SDMAX = 1000 ;
  size_t iters = 0 , i ;

  // allocate the temporary fitfunction for computing new steps 
  // down descent direction in the line search
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , 
					    Fit -> f.N ) ;

  // allocate the gradient
  double *grad = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  double *y = NULL , *t , alpha[ Fit -> Nlogic ] ;
  size_t Nsum = Fit -> f.N ;
  switch( Fit -> f.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = Fit -> f.N * Fit -> f.N ; break ;
  }
  y = malloc( Nsum * sizeof( double ) ) ;
  
  // get priors
  Fit -> f.Prior = Fit -> Prior ;
  f2.Prior = Fit -> Prior ;

  // evaluate the function, its first and second derivatives
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

  // set array of alphas to the fit parameters as an initial guess
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    alpha[i] = Fit -> f.fparams[i] ;
  }
  
  double chisq_diff = 10 ;
  while( chisq_diff > TOL && iters < SDMAX ) {

    // compute the descent direction ( - the gradient for SD! )
    for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
      size_t j , k ;
      t = y ;
      switch( Fit -> f.CORRFIT ) {
      case UNWEIGHTED :
	for( j = 0 ; j < Fit -> f.N ; j++ ) {
	  *t = -Fit -> f.df[i][j] * Fit -> f.f[j] ; t++ ;
	}
	break ;
      case UNCORRELATED :
	for( j = 0 ; j < Fit -> f.N ; j++ ) {
	  *t = -Fit -> f.df[i][j] * W[0][j] * Fit -> f.f[j] ; t++ ;
	}
	break ;
      case CORRELATED :
	for( j = 0 ; j < Fit -> f.N ; j++ ) {
	  for( k = 0 ; k < Fit -> f.N ; k++ ) {
	    *t = -Fit -> f.df[i][j] * W[j][k] * Fit -> f.f[k] ; t++ ;
	  }
	}
	break ;
      }
      grad[i] = kahan_summation( y , Nsum ) ;
      if( Fit -> f.Prior[i].Initialised == true ) {
	grad[i] -= ( Fit -> f.fparams[i] - Fit -> f.Prior[i].Val ) / 
	  ( Fit -> f.Prior[i].Err * Fit -> f.Prior[i].Err ) ;
      }
    }

    // do the line search, it is a good idea to make the search proportional
    // to the fparams[i] so that if there is a large difference between parameters
    // this gets taken into account
    for( i = 0 ; i < Fit -> Nlogic ; i++ ) {

      // line search best alpha
      alpha[i] = 1 ;
      alpha[i] = line_search( &f2 , Fit -> f , grad , grad ,
			      *Fit , data , W , i ,
			      alpha[i] ) ;

      #ifdef VERBOSE
      printf( "[SD] fparam :: %e \n" , Fit->f.fparams[i] ) ;
      printf( "[SD] line ap :: %e || grad %e \n" , alpha[i] , grad[i] ) ;
      #endif
      
      Fit -> f.fparams[i] += alpha[i] * grad[i] ;
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

  // free the temporary array
  if( y != NULL ) {
    free( y ) ;
  }

  // free the gradient
  free( grad ) ;

  // free the fitfunction
  free_ffunction( &f2 , f2.NPARAMS ) ;

  return iters ;
}

#undef BIG_GUESS
