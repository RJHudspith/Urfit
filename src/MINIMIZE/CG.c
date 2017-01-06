/**
   conjugate gradient routine
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h"

// verbose output
//#define VERBOSE

// depending on your function you might want to change these
#define BIG_GUESS (50) // the biggest guess for the backtracker

// cg iteration
int
cg_iter( struct fit_descriptor *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // iterations and CG maximum iterations
  size_t iters = 0 , i , j , k ;
  const size_t CGMAX = 5000 ;

  // allocate the fitfunction
  struct ffunction f2 = allocate_ffunction( fdesc -> Nlogic , 
					    fdesc -> f.N ) ;

  // set the guesses
  fdesc -> guesses( fdesc -> f.fparams , fdesc -> Nlogic ) ;

  // get priors
  fdesc -> set_priors( fdesc -> f.prior , fdesc -> f.err_prior ) ;

  // evaluate the function, its first and second derivatives
  fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ;
  fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;
  fdesc -> d2F( fdesc -> f.d2f , data , fdesc -> f.fparams ) ;
  fdesc -> f.chisq = compute_chisq( fdesc -> f , W , fdesc -> f.CORRFIT ) ;

  // allocate conjugate directions
  double *s      = malloc( fdesc -> Nlogic * sizeof( double* ) ) ;
  double *old_df = malloc( fdesc -> Nlogic * sizeof( double* ) ) ;

  // step down the gradient initially
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    old_df[i] = 0.0 ;
    // set derivatives
    for( j = 0 ; j < fdesc -> f.N ; j++ ) {
      switch( fdesc -> f.CORRFIT ) {
      case UNWEIGHTED :
	old_df[i] += -fdesc -> f.df[i][j] * fdesc -> f.f[j] ;
	break ;
      case UNCORRELATED :
	old_df[i] += -fdesc -> f.df[i][j] * W[j][j] * fdesc -> f.f[j] ;
	break ;
      case CORRELATED :
	for( k = 0 ; k < fdesc -> f.N ; k++ ) {
	  old_df[i] += -fdesc -> f.df[i][j] * W[j][k] * fdesc -> f.f[k] ;
	}
	break ;
      }
    }
    // subtract the prior if there is one
    if( fdesc -> f.prior[i] != UNINIT_FLAG ) {
      old_df[i] -= ( fdesc -> f.fparams[i] - fdesc -> f.prior[i] ) / 
	( fdesc -> f.err_prior[i] * fdesc -> f.err_prior[i] ) ;
    }
    s[i] = old_df[i] ;   // accumulate gradient sum
  }

  // line search the SD step
  double alpha = 0.0 ;
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    alpha = line_search( &f2 , fdesc -> f , old_df , s , 
			 *fdesc , data , W , i , BIG_GUESS ) ;
    fdesc -> f.fparams[i] += alpha * s[i] ;
  }

  double chisq_diff = 10 ;
  while( chisq_diff > TOL && iters < CGMAX ) {

    // initially compute the chisq
    const double chi = compute_chisq( fdesc -> f , W , 
				      fdesc -> f.CORRFIT ) ;
    
    // update f which has the gradient direction in
    fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ;
    fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;
    fdesc -> f.chisq = compute_chisq( fdesc -> f , W , fdesc -> f.CORRFIT ) ;
    chisq_diff = fabs( chi - fdesc -> f.chisq ) ;

    // compute beta using polyak - ribiere
    register double num = 0.0 , denom = 0.0 ;
    for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
      double newdf = 0.0 ;
      for( j = 0 ; j < fdesc -> f.N ; j++ ) {
	// switch the types of fit we are doing
	switch( fdesc -> f.CORRFIT ) {
	case UNWEIGHTED :
	  newdf += -fdesc -> f.df[i][j] * fdesc -> f.f[j] ;
	  break ;
	case UNCORRELATED : 
	  newdf += -fdesc -> f.df[i][j] * W[j][j] * fdesc -> f.f[j] ;
	  break ;
	case CORRELATED :
	  for( k = 0 ; k < fdesc -> f.N ; k++ ) {
	    newdf += -fdesc -> f.df[i][j] * W[j][k] * fdesc -> f.f[k] ;
	  }
	  break ;
	}
      }
      if( fdesc -> f.prior[i] != UNINIT_FLAG ) {
	newdf -= ( fdesc -> f.fparams[i] - fdesc -> f.prior[i] ) / 
	  ( fdesc -> f.err_prior[i] * fdesc -> f.err_prior[i] ) ;
      }
      num   += newdf * ( newdf - old_df[i] ) ;
      denom += old_df[i] * old_df[i] ;
      old_df[i] = newdf ; // reset the old direction to be the new one
    }
    const double beta = fmax( 0 , num / denom ) ; 

    #ifdef VERBOSE
    printf( "[CG] BETA :: %e ( %e / %e ) \n" , beta , num , denom ) ;
    #endif

    // update conjugate directions "s"
    for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
      s[i] = old_df[i] + beta * s[i] ;
    }

    for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
      
      // perform a backtracking line search
      alpha = line_search( &f2 , fdesc -> f , old_df , s , 
			   *fdesc , data , W , i , BIG_GUESS ) ;

      #ifdef VERBOSE
      printf( "[CG] NEW ALPHA = %e \n" , alpha ) ;
      #endif

      // set x = x + \alpha * s
      fdesc -> f.fparams[i] += alpha * s[i] ;

      #ifdef VERBOSE
      printf( "[CG] NEPARAMS :: %f \n" , fdesc -> f.fparams[i] ) ;
      #endif
    }
    iters++ ;
  }

  // tell us how many iterations we hit
  if( iters == CGMAX ) {
    printf( "\n[CG] stopped by max iterations %zu \n" , iters ) ;
  } else {
    printf( "\n[CG] FINISHED in %zu iterations \n" , iters ) ;
  }

  
  printf( "[CG] chisq :: %e | Diff -> %e \n\n" , fdesc -> f.chisq , chisq_diff ) ;
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    printf( "PARAMS :: %f \n" , fdesc -> f.fparams[i] ) ;
  }

  // free the directions
  free( s ) ;
  free( old_df ) ;

  return iters ;
}

// clean up the defines
#ifdef VERBOSE
  #undef VERBOSE
#endif
#undef BIG_GUESS

