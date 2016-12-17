/**
   @file LM.c
   @brief levenberg-marquardt algorithm

   Marquardt - levenberg algorithm
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"

//#define VERBOSE

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
	bp += f.df[p][i] * W[i][i] * f.f[i] ;
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
    if( f.prior[p] != UNINIT_FLAG ) {
      bp += ( f.fparams[p] - f.prior[p] ) / 
	( f.err_prior[p] * f.err_prior[p] ) ;
    }

    #ifdef VERBOSE
    printf( "[ML] beta[%zu] %f \n" , p , bp ) ;
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
	    + f.d2f[q+f.NPARAMS*p][i] * f.f[i] 
	    #endif
	    ;
	}
	break ;
      case UNCORRELATED :
	for( i = 0 ; i < f.N ; i++ ) {
	  apq += 
	    W[i][i] * ( f.df[p][i] * f.df[q][i] 
                        #ifdef WITH_D2_DERIVS
			+ f.d2f[q+f.NPARAMS*p][i] * f.f[i]
			#endif
			)
	    ;
	}
	break ;
      case CORRELATED :
	for( i = 0 ; i < f.N ; i++ ) {
	  for( j = 0 ; j < f.N ; j++ ) {
	    apq += 
	      W[i][j] * ( f.df[p][i] * f.df[q][j] 
			  #ifdef WITH_D2_DERIVS
			  + f.d2f[q+f.NPARAMS*p][i] * f.f[j] 
			  #endif
			  ) 
	      ;
	  }
	}
	break ;
      }

      // second derivatives acting on the prior
      if( p == q ) {
	if( f.prior[p] != UNINIT_FLAG ) {
	  apq += 1.0 / 
	    ( f.err_prior[p] * f.err_prior[p] ) ;
	}
      }

      #ifdef VERBOSE
      printf( "[ML] alpha[%zu,%zu] %f \n" , p , q  , -apq ) ;
      #endif

      gsl_matrix_set( alpha , p , q , -apq ) ;
    }
    // end of setup
  }
  return GSL_SUCCESS ;
}

// perform trial - ML step return chisq
double
ml_step( struct ffunction *f , 
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
    printf( "[ML] LU decomp broken?\n" ) ;
    return !GSL_SUCCESS ;
  }
  if( gsl_linalg_LU_solve( alpha_new , perm , beta , delta ) != GSL_SUCCESS ) {
    printf( "[ML] LU solve broken?\n" ) ;
    return !GSL_SUCCESS ;
  }
  // update fitparams
  for( i = 0 ; i < f -> NPARAMS ; i++ ) {
    f -> fparams[i] = old_params[i] + gsl_vector_get( delta , i ) ;
    #ifdef VERBOSE
    printf( "[ML] NEW PARAMS :: %f \n" , f -> fparams[i] ) ;
    #endif
  }
  // evaluate these functions
  fdesc.F( f -> f , data , f -> fparams ) ;
  // compute new chisq
  return compute_chisq( *f , W , f -> CORRFIT ) ;
}

// perform marquardt - levenberg updates
int
ml_iter( struct fit_descriptor *fdesc ,
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
  fdesc -> guesses( fdesc -> f.fparams ) ;

  // get priors
  fdesc -> set_priors( fdesc -> f.prior , fdesc -> f.err_prior ) ;

  // allocate alpha, beta, delta and permutation matrices
  gsl_matrix *alpha = gsl_matrix_alloc( fdesc -> NPARAMS , 
					fdesc -> NPARAMS ) ;
  gsl_matrix *alpha_new = gsl_matrix_alloc( fdesc -> NPARAMS , 
					    fdesc -> NPARAMS ) ;
  gsl_vector *beta = gsl_vector_alloc( fdesc -> NPARAMS ) ;
  gsl_vector *delta = gsl_vector_alloc( fdesc -> NPARAMS ) ;
  gsl_permutation *perm = gsl_permutation_alloc( fdesc -> NPARAMS ) ;

  // set the old parameters
  double old_params[ fdesc -> NPARAMS ] ;
  for( i = 0 ; i < fdesc -> NPARAMS ; i++ ) {
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

    // TODO :: line search or trust region for Lambda!

    double new_chisq = ml_step( &fdesc -> f , old_params , alpha_new , delta , 
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
      for( i = 0 ; i < fdesc -> NPARAMS ; i++ ) {
	old_params[ i ] = fdesc -> f.fparams[ i ] ;
      }
      // update derivatives and alpha and beta
      fdesc -> dF( fdesc -> f.df , data , fdesc -> f.fparams ) ;
      #ifdef WITH_D2_DERIVS
      fdesc -> d2F( fdesc -> f.d2f , data , fdesc -> f.fparams ) ;
      #endif
      get_alpha_beta( alpha , beta , fdesc -> f , W ) ;
    } else {
      for( i = 0 ; i < fdesc -> NPARAMS ; i++ ) {
	fdesc -> f.fparams[ i ] = old_params[ i ] ;
      }
      fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ; // reset f
      Lambda *= fac ;
    }

    // increment the number of iterations
    iters++ ;
  }

  // tell us how many iterations we hit
  if( iters == LMMAX ) {
    printf( "\n[ML] stopped by max iterations %zu \n" , iters ) ;
  }

#ifdef VERBOSE
  printf( "\n[ML] FINISHED in %zu iterations \n" , iters ) ;
  printf( "[ML] chisq :: %e \n\n" , fdesc -> f.chisq ) ;
  // tell us the fit parameters
  for( i = 0 ; i < fdesc -> NPARAMS ; i++ ) {
    printf( "PARAMS :: %f \n" , fdesc -> f.fparams[i] ) ;
  }
#endif

  // free gsl stuff
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( alpha_new ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  return iters ;
}
