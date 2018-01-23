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
//#define WITH_D2_DERIVS

//#define VERBOSE
#define TRUST_REGION

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
	 double *pred ,
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

  // compute the predicted chi^2
  register double loc_sum = f -> chisq ;
  for( i = 0 ; i < f -> NPARAMS ; i++ ) {
    loc_sum += gsl_vector_get( delta , i ) * gsl_vector_get( beta , i ) ;
    size_t j ;
    for( j = 0 ; j < f -> NPARAMS ; j++ ) {
      loc_sum +=
	0.5 * gsl_vector_get( delta , i ) *
	gsl_matrix_get( alpha , i , j ) *
	gsl_vector_get( delta , j ) ;
    }
  }
  *pred = loc_sum ;
  
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
  const size_t LMMAX = 1000 ; 

  // lambda multiplying factor
  const double fac = 5. ;

  double chisq_diff = 10 , Lambda = 0.5 ;
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

#ifdef TRUST_REGION
  const double Lambda_max = 1E10 ;
#endif

  // loop until chisq evens out
  double pred = 1.0 ;
  while( chisq_diff > TOL && iters < LMMAX ) {

    const double new_chisq = lm_step( &Fit -> f , &pred , old_params ,
				      alpha_new , delta , 
				      perm , alpha , beta , *Fit , 
				      data , W , Lambda ) ;

    const double ratio = ( Fit -> f.chisq - new_chisq ) /
      ( Fit -> f.chisq - pred ) ;

#ifdef TRUST_REGION

    const double eta1 = 0.25 , eta2 = 1.0 ;
    register double frob = 0.0 ;
    for( i = 0 ; i < Fit -> f.NPARAMS ; i++ ) {
      frob += gsl_vector_get( beta , i ) * gsl_vector_get( beta , i ) ;
    }
    frob = sqrt( frob ) ;

    #ifdef VERBOSE
    printf( "Ratio %e | Lambda_prev %e \n" , ratio , Lambda ) ;
    printf( "%zu Frob %e \n" , iters , frob ) ;
    #endif

    // update trust region radius
    if( ratio < eta1 ) {
      Lambda = fmin( Lambda * 4. , Lambda_max ) ;
    } else if( ratio > eta2 ) {
      Lambda = Lambda / 4. ;
    }
    
    // update as we are in the trust region
    if( ratio > eta1 ) {
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
    }
    #ifdef VERBOSE
    printf( "Chidiff %e \n" , chisq_diff ) ;
    printf( "Lambda %e \n" , Lambda ) ;
    #endif

    // if the trust region keeps expanding then we just quit
    if( Lambda >= Lambda_max ) {
      //fprintf( stderr , "Lambda == Lambda_max\n" ) ;
      break ;
    }
#else // resort to original version
    
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
#endif

    // if lambda misses a minimum tell us
    if( Lambda < 1E-32 || Lambda > 1E32 ) {
      printf( "[LM] Lambda is out of bounds %e \n" , Lambda ) ;
      iters = LMMAX ;
      break ;
    }
    
    // increment the number of iterations
    iters++ ;
  }
  
  #ifdef VERBOSE
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
  #endif

  // free gsl stuff
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( alpha_new ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  return iters ;
}
