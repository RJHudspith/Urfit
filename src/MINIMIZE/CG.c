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
cg_iter( void *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // point to the fit descriptor
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  
  // iterations and CG maximum iterations
  size_t iters = 0 , i , j , k ;
  const size_t CGMAX = 5000 ;

  // allocate the fitfunction
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , 
					    Fit -> f.N ) ;

  // get priors
  Fit -> f.Prior = Fit -> Prior ;
  f2.Prior = Fit -> Prior ;

  // evaluate the function, its first and second derivatives
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
  //Fit -> d2F( Fit -> f.d2f , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

  // allocate conjugate directions
  double *s      = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  double *old_df = malloc( Fit -> Nlogic * sizeof( double ) ) ;

  double *y = NULL , *t ;
  size_t Nsum = Fit -> f.N ;
  switch( Fit -> f.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = Fit -> f.N * Fit -> f.N ; break ;
  }
  y = malloc( Nsum * sizeof( double ) ) ;

  // step down the gradient initially
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    t = y ;
    // set derivatives
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
    old_df[i] = kahan_summation( y , Nsum ) ;
    // subtract the prior if there is one
    if( Fit -> f.Prior[i].Initialised == true ) {
      old_df[i] -= ( Fit -> f.fparams[i] - Fit -> f.Prior[i].Val ) / 
	( Fit -> f.Prior[i].Err * Fit -> f.Prior[i].Err ) ;
    }
    s[i] = old_df[i] ;   // accumulate gradient sum
  }

  // line search the SD step
  double alpha = 0.0 ;
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    // line search best alpha
    alpha = line_search( &f2 , Fit -> f , old_df , s ,
			 *Fit , data , W , i ,
			 Fit -> f.fparams[i] ) ;
    Fit -> f.fparams[i] += alpha * s[i] ;
  }

  double chisq_diff = 10 ;
  while( chisq_diff > TOL && iters < CGMAX ) {

    // initially compute the chisq
    const double chi = compute_chisq( Fit -> f , W , 
				      Fit -> f.CORRFIT ) ;
    
    // update f which has the gradient direction in
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
    Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
    Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
    chisq_diff = fabs( chi - Fit -> f.chisq ) ;

    // compute beta using polyak - ribiere
    register double num = 0.0 , denom = 0.0 ;
    for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
      t = y ;
      // switch the types of fit we are doing
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
      double newdf = kahan_summation( y , Nsum ) ;
      
      // subtract the prior if there is one
      if( Fit -> f.Prior[i].Initialised == true ) {
        newdf -= ( Fit -> f.fparams[i] - Fit -> f.Prior[i].Val ) / 
	  ( Fit -> f.Prior[i].Err * Fit -> f.Prior[i].Err ) ;
      }
      num   += newdf * ( newdf - old_df[i] ) ;
      denom += old_df[i] * old_df[i] ;
      old_df[i] = newdf ; // reset the old direction to be the new one
    }
    const double beta = fmax( 0 , num / denom ) ;

    // update conjugate directions "s"
    for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
      s[i] = old_df[i] + beta * s[i] ;
    }

    // perform a backtracking line search
    for( i = 0 ; i < Fit -> Nlogic ; i++ ) {

      alpha = line_search( &f2 , Fit -> f , old_df , s ,
			   *Fit , data , W , i ,
			   Fit -> f.fparams[i] ) ;

      #ifdef VERBOSE
      printf( "[CG] NEW ALPHA = %e \n" , alpha ) ;
      #endif

      // set x = x + \alpha * s
      Fit -> f.fparams[i] += alpha * s[i] ;

      #ifdef VERBOSE
      printf( "[CG] NEPARAMS :: %f \n" , Fit -> f.fparams[i] ) ;
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

  printf( "[CG] chisq :: %e | Diff -> %e \n\n" , Fit -> f.chisq , chisq_diff ) ;
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    printf( "PARAMS :: %f \n" , Fit -> f.fparams[i] ) ;
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

