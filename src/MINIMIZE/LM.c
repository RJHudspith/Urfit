/**
   @file LM.c
   @brief levenberg-marquardt algorithm

   Marquardt - levenberg algorithm
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "summation.h"

// use second derivs? NRC says that they can be harmful
// not using them makes the code run faster sometimes
#define WITH_D2_DERIVS 

// get the matrices alpha and beta
static int
get_alpha_beta( gsl_matrix *alpha ,
		gsl_vector *beta ,
		const struct ffunction f ,
		const double **W )
{
  double *y = NULL , *t ;
  size_t p , q = 0 , i , j , Nsum = f.N ;

  // allocate the y-data
  switch( f.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = f.N * f.N ; break ;
  }
  y = malloc( Nsum * sizeof( double ) ) ;
  
  // compute beta
  for( p = 0 ; p < f.NPARAMS ; p++ ) {

    // compute beta
    t = y ;
    switch( f.CORRFIT ) {
    case UNWEIGHTED : 
      for( i = 0 ; i < f.N ; i++ ) {
	*t = f.df[p][i] * f.f[i] ; t++ ;
      }
      break ;
    case UNCORRELATED :
      for( i = 0 ; i < f.N ; i++ ) {
	*t = f.df[p][i] * W[0][i] * f.f[i] ; t++ ;
      }
      break ;
    case CORRELATED :
      for( i = 0 ; i < f.N ; i++ ) {
	for( j = 0 ; j < f.N ; j++ ) {
	  *t = f.df[p][i] * W[i][j]  * f.f[j] ; t++ ;
	}
      }
      break ;
    }
    register double bp = kahan_summation( y , Nsum ) ;

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

      // compute alpha
      t = y ;
      switch( f.CORRFIT ) {
      case UNWEIGHTED : 
	for( i = 0 ; i < f.N ; i++ ) {
	  *t = 
	    f.df[p][i] * f.df[q][i] 
	    #ifdef WITH_D2_DERIVS
	    + f.d2f[q+f.NPARAMS*p][i] * f.f[i] 
	    #endif
	    ; t++ ;
	}
	break ;
      case UNCORRELATED :
	for( i = 0 ; i < f.N ; i++ ) {
	  *t = 
	    W[0][i] * ( f.df[p][i] * f.df[q][i] 
                        #ifdef WITH_D2_DERIVS
			+ f.d2f[q+f.NPARAMS*p][i] * f.f[i]
			#endif
			) ; t++ ;
	}
	break ;
      case CORRELATED :
	for( i = 0 ; i < f.N ; i++ ) {
	  for( j = 0 ; j < f.N ; j++ ) {
	    *t = 
	      W[i][j] * ( f.df[p][i] * f.df[q][j] 
			  #ifdef WITH_D2_DERIVS
			  + f.d2f[q+f.NPARAMS*p][i] * f.f[j] 
			  #endif
			  ) ; t++ ;
	  }
	}
	break ;
      }
      register double apq = kahan_summation( y , Nsum ) ;
      
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

  // clean up the temporary data
  if( y != NULL ) {
    free( y ) ;
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
lm_iter( void *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // point to the fit function
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  
  // set maximum iterations
  const size_t LMMAX = 10000 ; 

  // lambda multiplying factor
  const double fac = 10 ;

  double chisq_diff = 10 , Lambda = 0.1 ;
  size_t iters = 0 , i ;

  // get priors
  Fit -> f.Prior = Fit -> Prior ;

  // allocate alpha, beta, delta and permutation matrices
  gsl_matrix *alpha     = gsl_matrix_alloc( Fit -> Nlogic , 
					    Fit -> Nlogic ) ;
  gsl_matrix *alpha_new = gsl_matrix_alloc( Fit -> Nlogic , 
					    Fit -> Nlogic ) ;
  gsl_vector *beta  = gsl_vector_alloc( Fit -> Nlogic ) ;
  gsl_vector *delta = gsl_vector_alloc( Fit -> Nlogic ) ;
  gsl_permutation *perm = gsl_permutation_alloc( Fit -> Nlogic ) ;

  // set the old parameters
  double old_params[ Fit -> Nlogic ] ;
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    old_params[ i ] = Fit -> f.fparams[ i ] ;
  }

  // evaluate the function, its first and second derivatives
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
  #ifdef WITH_D2_DERIVS
  Fit -> d2F( Fit -> f.d2f , data , Fit -> f.fparams ) ;
  #endif
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

  // get alpha and beta for the new set of f
  get_alpha_beta( alpha , beta , Fit -> f , W ) ;

  // loop until chisq evens out
  while( chisq_diff > TOL ) {

    double new_chisq = lm_step( &Fit -> f , old_params , alpha_new , delta , 
				perm , alpha , beta , *Fit , 
				data , W , Lambda ) ;    

    #ifdef VERBOSE
    printf( "[ML] chis :: %f %f %e \n" , new_chisq , 
	    Fit -> f.chisq , fabs( Fit -> f.chisq - new_chisq ) ) ;
    #endif

    // update lambda shrinking the size if we are close to a solution
    if( new_chisq <= Fit -> f.chisq && iters < LMMAX ) {
      // update lambda
      Lambda /= fac ;
      chisq_diff = fabs( Fit -> f.chisq - new_chisq ) ;
      Fit -> f.chisq = new_chisq ;
      for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
	old_params[ i ] = Fit -> f.fparams[ i ] ;
      }
      // update derivatives and alpha and beta
      Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
      #ifdef WITH_D2_DERIVS
      Fit -> d2F( Fit -> f.d2f , data , Fit -> f.fparams ) ;
      #endif
      get_alpha_beta( alpha , beta , Fit -> f , W ) ;
    } else {
      for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
	Fit -> f.fparams[ i ] = old_params[ i ] ;
      }
      Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ; // reset f
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
  
  //#ifdef VERBOSE
  // tell us how many iterations we hit
  if( iters == LMMAX ) {
    printf( "\n[LM] stopped by max iterations %zu -> Chidiff %e\n" , iters , chisq_diff ) ;
  } else {
    printf( "\n[LM] FINISHED in %zu iterations \n" , iters ) ;
  }
  
  printf( "[LM] chisq :: %e \n\n" , Fit -> f.chisq ) ;
  // tell us the fit parameters
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    printf( "PARAMS :: %f \n" , Fit -> f.fparams[i] ) ;
  }
  //#endif

  // free gsl stuff
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( alpha_new ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  return iters ;
}
