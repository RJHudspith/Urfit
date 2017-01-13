/**
   @file LM.c
   @brief levenberg-marquardt algorithm

   Marquardt - levenberg algorithm
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"

// use second derivs? NRC says that they can be harmful
// not using them makes the code run faster sometimes
//#define WITH_D2_DERIVS 

// get the matrices alpha and beta
static int
get_alpha_beta( gsl_matrix *alpha ,
		gsl_vector *beta ,
		const struct ffunction f ,
		const double **W )
{
  // compute beta
  size_t p , q = 0 , i , j ;
  for( p = 0 ; p < f.NPARAMS ; p++ ) {

    // compute beta
    register double bp = 0.0 ;

    switch( f.CORRFIT ) {
    case UNWEIGHTED : 
      for( i = 0 ; i < f.N ; i++ ) {
	bp += f.df[p][i] * f.f[i] ;
      }
      break ;
    case UNCORRELATED :
      for( i = 0 ; i < f.N ; i++ ) {
	bp += f.df[p][i] * W[0][i] * f.f[i] ;
      }
      break ;
    case CORRELATED :
      for( i = 0 ; i < f.N ; i++ ) {
	for( j = 0 ; j < f.N ; j++ ) {
	  bp += f.df[p][i] * W[i][j]  * f.f[j] ;
	}
      }
      break ;
    }

    // add the priors if they have been set
    if( f.Prior[p].Initialised == true ) {
      bp += ( f.fparams[p] - f.Prior[p].Val ) / 
	( f.Prior[p].Err * f.Prior[p].Err ) ;
    }

    #ifdef VERBOSE
    printf( "[LM] beta[%zu] %f \n" , p , bp ) ;
    #endif

    // set beta[p]
    gsl_vector_set( beta , p , bp ) ;

    // compute alpha[p][q].loop q only need to do top half
    // as the matrix "alpha" is symmetric
    for( q = p ; q < f.NPARAMS ; q++ ) {

      // set alpha
      register double apq = 0.0 ;

      // compute alpha
      switch( f.CORRFIT ) {
      case UNWEIGHTED : 
	for( i = 0 ; i < f.N ; i++ ) {
	  apq += 
	    f.df[p][i] * f.df[q][i] 
	    #ifdef WITH_D2_DERIVS
	    + f.d2f[q+f.Nlogic*p][i] * f.f[i] 
	    #endif
	    ;
	}
	break ;
      case UNCORRELATED :
	for( i = 0 ; i < f.N ; i++ ) {
	  apq += 
	    W[0][i] * ( f.df[p][i] * f.df[q][i] 
                        #ifdef WITH_D2_DERIVS
			+ f.d2f[q+f.Nlogic*p][i] * f.f[i]
			#endif
			) ;
	}
	break ;
      case CORRELATED :
	for( i = 0 ; i < f.N ; i++ ) {
	  for( j = 0 ; j < f.N ; j++ ) {
	    apq += 
	      W[i][j] * ( f.df[p][i] * f.df[q][j] 
			  #ifdef WITH_D2_DERIVS
			  + f.d2f[q+f.Nlogic*p][i] * f.f[j] 
			  #endif
			  ) ;
	  }
	}
	break ;
      }

      // second derivatives acting on the prior
      if( p == q ) {
	if( f.Prior[p].Initialised == true ) {
	  apq += 1.0 / 
	    ( f.Prior[p].Err * f.Prior[p].Err ) ;
	}
      }

      #ifdef VERBOSE
      printf( "[LM] alpha[%zu,%zu] %f \n" , p , q  , -apq ) ;
      #endif

      gsl_matrix_set( alpha , p , q , -apq ) ;
    }
    // end of setup
  }
  return GSL_SUCCESS ;
}

// perform trial - ML step return chisq
double
lm_step( struct ffunction *f , 
	 double *old_params ,
	 gsl_matrix *alpha_new ,
	 gsl_vector *delta ,
	 gsl_permutation *perm ,
	 const gsl_matrix *alpha ,
	 const gsl_vector *beta ,
	 const struct fit_descriptor fdesc ,
	 const void *data ,
	 const double **W ,
	 const double Lambda )
{
  int signum ;
  size_t i , j ;
  // compute a new alpha by adjusting the diagonal
  for( i = 0 ; i < f -> NPARAMS ; i++ ) {
    for( j = i+1 ; j < f -> NPARAMS ; j++ ) {
      const double ref = gsl_matrix_get( alpha , i , j ) ;
      gsl_matrix_set( alpha_new , i , j , ref ) ;
      gsl_matrix_set( alpha_new , j , i , ref ) ;
    }
    gsl_matrix_set( alpha_new , i , i , 
		    ( 1.0 + Lambda ) * gsl_matrix_get( alpha , i , i ) ) ;
  }
  // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
  if( gsl_linalg_LU_decomp( alpha_new , perm , &signum ) != GSL_SUCCESS ) {
    printf( "[LM] LU decomp broken?\n" ) ;
    return !GSL_SUCCESS ;
  }
  if( gsl_linalg_LU_solve( alpha_new , perm , beta , delta ) != GSL_SUCCESS ) {
    printf( "[LM] LU solve broken?\n" ) ;
    return !GSL_SUCCESS ;
  }
  // update fitparams
  for( i = 0 ; i < f -> NPARAMS ; i++ ) {
    f -> fparams[i] = old_params[i] + gsl_vector_get( delta , i ) ;
    #ifdef VERBOSE
    printf( "[LM] NEW PARAMS :: %f \n" , f -> fparams[i] ) ;
    #endif
  }
  // evaluate these functions
  fdesc.F( f -> f , data , f -> fparams ) ;
  // compute new chisq
  return compute_chisq( *f , W , f -> CORRFIT ) ;
}

// perform marquardt - levenberg updates
int
lm_iter( struct fit_descriptor *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // allocate the fitfunction
  const size_t LMMAX = 10000 ; 

  // lambda multiplying factor
  const double fac = 10 ;

  double chisq_diff = 10 , Lambda = 0.1 ;
  size_t iters = 0 , i ;

  // some guesses
  fdesc -> guesses( fdesc -> f.fparams , fdesc -> Nlogic ) ;

  // get priors
  fdesc -> f.Prior = fdesc -> Prior ;

  // allocate alpha, beta, delta and permutation matrices
  gsl_matrix *alpha     = gsl_matrix_alloc( fdesc -> Nlogic , 
					    fdesc -> Nlogic ) ;
  gsl_matrix *alpha_new = gsl_matrix_alloc( fdesc -> Nlogic , 
					    fdesc -> Nlogic ) ;
  gsl_vector *beta  = gsl_vector_alloc( fdesc -> Nlogic ) ;
  gsl_vector *delta = gsl_vector_alloc( fdesc -> Nlogic ) ;
  gsl_permutation *perm = gsl_permutation_alloc( fdesc -> Nlogic ) ;

  // set the old parameters
  double old_params[ fdesc -> Nlogic ] ;
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    old_params[ i ] = fdesc -> f.fparams[ i ] ;
  }

  // evaluate the function, its first and second derivatives
  fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ;
  fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;
  #ifdef WITH_D2_DERIVS
  fdesc -> d2F( fdesc -> f.d2f , data , fdesc -> f.fparams ) ;
  #endif
  fdesc -> f.chisq = compute_chisq( fdesc -> f , W , fdesc -> f.CORRFIT ) ;

  // get alpha and beta for the new set of f
  get_alpha_beta( alpha , beta , fdesc -> f , W ) ;

  // loop until chisq evens out
  while( chisq_diff > TOL ) {

    double new_chisq = lm_step( &fdesc -> f , old_params , alpha_new , delta , 
				perm , alpha , beta , *fdesc , 
				data , W , Lambda ) ;    

    #ifdef VERBOSE
    printf( "[ML] chis :: %f %f %e \n" , new_chisq , 
	    fdesc -> f.chisq , fabs( fdesc -> f.chisq - new_chisq ) ) ;
    #endif

    // update lambda shrinking the size if we are close to a solution
    if( new_chisq <= fdesc -> f.chisq && iters < LMMAX ) {
      // update lambda
      Lambda /= fac ;
      chisq_diff = fabs( fdesc -> f.chisq - new_chisq ) ;
      fdesc -> f.chisq = new_chisq ;
      for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
	old_params[ i ] = fdesc -> f.fparams[ i ] ;
      }
      // update derivatives and alpha and beta
      fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;
      #ifdef WITH_D2_DERIVS
      fdesc -> d2F( fdesc -> f.d2f , data , fdesc -> f.fparams ) ;
      #endif
      get_alpha_beta( alpha , beta , fdesc -> f , W ) ;
    } else {
      for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
	fdesc -> f.fparams[ i ] = old_params[ i ] ;
      }
      fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ; // reset f
      Lambda *= fac ;
    }

    // if lambda misses a minimum tell us
    if( Lambda < 1E-45 || Lambda > 1E45 ) {
      printf( "[LM] Lambda is out of bounds %e \n" , Lambda ) ;
      iters = LMMAX ;
      break ;
    }
    
    // increment the number of iterations
    iters++ ;
  }

  // tell us how many iterations we hit
  if( iters == LMMAX ) {
    printf( "\n[LM] stopped by max iterations %zu -> Chidiff %e\n" , iters , chisq_diff ) ;
  } else {
    printf( "\n[LM] FINISHED in %zu iterations \n" , iters ) ;
  }
  
  printf( "[LM] chisq :: %e \n\n" , fdesc -> f.chisq ) ;
  // tell us the fit parameters
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    printf( "PARAMS :: %f \n" , fdesc -> f.fparams[i] ) ;
  }

  // free gsl stuff
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( alpha_new ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  return iters ;
}
