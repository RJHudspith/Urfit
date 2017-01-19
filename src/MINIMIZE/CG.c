/**
   conjugate gradient routine
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h"
#include "summation.h"

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
  fdesc -> f.Prior = fdesc -> Prior ;
  f2.Prior = fdesc -> Prior ;

  // evaluate the function, its first and second derivatives
  fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ;
  fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;
  fdesc -> d2F( fdesc -> f.d2f , data , fdesc -> f.fparams ) ;
  fdesc -> f.chisq = compute_chisq( fdesc -> f , W , fdesc -> f.CORRFIT ) ;

  // allocate conjugate directions
  double *s      = malloc( fdesc -> Nlogic * sizeof( double ) ) ;
  double *old_df = malloc( fdesc -> Nlogic * sizeof( double ) ) ;

  double *y = NULL , *t ;
  size_t Nsum = fdesc -> f.N ;
  switch( fdesc -> f.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = fdesc -> f.N * fdesc -> f.N ; break ;
  }
  y = malloc( Nsum * sizeof( double ) ) ;

  // step down the gradient initially
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    t = y ;
    // set derivatives
    switch( fdesc -> f.CORRFIT ) {
    case UNWEIGHTED :
      for( j = 0 ; j < fdesc -> f.N ; j++ ) {
	*t = -fdesc -> f.df[i][j] * fdesc -> f.f[j] ; t++ ;
      }
      break ;
    case UNCORRELATED :
      for( j = 0 ; j < fdesc -> f.N ; j++ ) {
        *t = -fdesc -> f.df[i][j] * W[0][j] * fdesc -> f.f[j] ; t++ ;
      }
      break ;
    case CORRELATED :
      for( j = 0 ; j < fdesc -> f.N ; j++ ) {
	for( k = 0 ; k < fdesc -> f.N ; k++ ) {
	  *t = -fdesc -> f.df[i][j] * W[j][k] * fdesc -> f.f[k] ; t++ ;
	}
      }
      break ;
    }
    old_df[i] = kahan_summation( y , Nsum ) ;
    // subtract the prior if there is one
    if( fdesc -> f.Prior[i].Initialised == true ) {
      old_df[i] -= ( fdesc -> f.fparams[i] - fdesc -> f.Prior[i].Val ) / 
	( fdesc -> f.Prior[i].Err * fdesc -> f.Prior[i].Err ) ;
    }
    s[i] = old_df[i] ;   // accumulate gradient sum
  }

  // line search the SD step
  double alpha = 0.0 ;
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    // line search best alpha
    alpha = line_search( &f2 , fdesc -> f , old_df , s ,
			 *fdesc , data , W , i ,
			 fdesc -> f.fparams[i] ) ;
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
      t = y ;
      // switch the types of fit we are doing
      switch( fdesc -> f.CORRFIT ) {
      case UNWEIGHTED :
	for( j = 0 ; j < fdesc -> f.N ; j++ ) {
	  *t = -fdesc -> f.df[i][j] * fdesc -> f.f[j] ; t++ ;
	}
	break ;
      case UNCORRELATED :
	for( j = 0 ; j < fdesc -> f.N ; j++ ) {
	  *t = -fdesc -> f.df[i][j] * W[0][j] * fdesc -> f.f[j] ; t++ ;
	}
	break ;
      case CORRELATED :
	for( j = 0 ; j < fdesc -> f.N ; j++ ) {
	  for( k = 0 ; k < fdesc -> f.N ; k++ ) {
	    *t = -fdesc -> f.df[i][j] * W[j][k] * fdesc -> f.f[k] ; t++ ;
	  }
	}
	break ;
      }
      double newdf = kahan_summation( y , Nsum ) ;
      
      // subtract the prior if there is one
      if( fdesc -> f.Prior[i].Initialised == true ) {
        newdf -= ( fdesc -> f.fparams[i] - fdesc -> f.Prior[i].Val ) / 
	  ( fdesc -> f.Prior[i].Err * fdesc -> f.Prior[i].Err ) ;
      }
      num   = newdf * ( newdf - old_df[i] ) ;
      denom = old_df[i] * old_df[i] ;
      old_df[i] = newdf ; // reset the old direction to be the new one

      const double beta = fmax( 0 , num / denom ) ;

      // update conjugate directions "s"
      s[i] = old_df[i] + beta * s[i] ;
    }

    // perform a backtracking line search
    for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {

      alpha = line_search( &f2 , fdesc -> f , old_df , s ,
			   *fdesc , data , W , i ,
			   fdesc -> f.fparams[i] ) ;

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
  if( s != NULL ) {
    free( s ) ;
  }

  if( old_df != NULL ) {
    free( old_df ) ;
  }

  if( y != NULL ) {
    free( y ) ;
  }

  return iters ;
}

// clean up the defines
#ifdef VERBOSE
  #undef VERBOSE
#endif
#undef BIG_GUESS

