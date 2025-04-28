/**
   @file LM.c
   @brief levenberg-marquardt algorithm

   Marquardt - levenberg algorithm
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "summation.h"
#include <gsl/gsl_errno.h>

// use second derivs? NRC says that they can be harmful
// not using them makes the code run faster sometimes
//#define WITH_D2_DERIVS

//#define VERBOSE
#define LMSVD

// allocate alpha, beta, delta and permutation matrices
// and other storage types for the LM
struct lmstep {
  gsl_matrix *alpha ;
  gsl_matrix *alpha_new ;
  gsl_vector *beta ; 
  gsl_vector *delta ; 
  gsl_permutation *perm ;
  double *old_params ;
  double *y ;
  double pred ;
  size_t Nsum ;
#ifdef LMSVD
  gsl_matrix *V ;
  gsl_vector *S ;
  gsl_vector *work ;
#endif
} ;

// get the matrices alpha and beta
static int
get_alpha_beta( struct lmstep *LM ,
		const struct ffunction f ,
		const double **W )
{
  double *t ;
  size_t p , q = 0 , i , j ;
  // compute beta gradient of \chi^2 function
  for( p = 0 ; p < f.NPARAMS ; p++ ) {
    t = LM -> y ;
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
    register double bp = kahan_summation( LM -> y , LM -> Nsum ) ;
    // add the priors if they have been set
    if( f.Prior[p].Initialised == true ) {
      bp += ( f.fparams[p] - f.Prior[p].Val ) / 
	( f.Prior[p].Err * f.Prior[p].Err ) ;
    }
    #ifdef VERBOSE
    printf( "[LM] beta[%zu] %e \n" , p , bp ) ;
    #endif
    // set beta[p]
    gsl_vector_set( LM -> beta , p , bp ) ;
    // compute alpha[p][q].loop q only need to do top half
    // as the matrix "alpha" is symmetric
    for( q = p ; q < f.NPARAMS ; q++ ) {
      // compute alpha
      t = LM -> y ;
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
      register double apq = kahan_summation( LM -> y , LM -> Nsum ) ;
      // second derivatives acting on the prior
      if( p == q ) {
	if( f.Prior[p].Initialised == true ) {
	  apq += 1.0 / 
	    ( f.Prior[p].Err * f.Prior[p].Err ) ;
	}
      }
      #ifdef VERBOSE
      printf( "[LM] alpha[%zu,%zu] %e \n" , p , q  , -apq ) ;
      #endif
      gsl_matrix_set( LM -> alpha , p , q , -apq ) ;
    }
  }
  
  return GSL_SUCCESS ;
}

// perform trial - ML step return chisq
static double
lm_step( struct ffunction *f ,
	 struct lmstep *LM ,
	 const struct fit_descriptor fdesc ,
	 const void *data ,
	 const double **W ,
	 const double Lambda )
{
  size_t i , j ;
  // compute a new alpha by adjusting the correctly-scaled diagonal
  for( i = 0 ; i < f -> NPARAMS ; i++ ) {
    for( j = i+1 ; j < f -> NPARAMS ; j++ ) {
      const double ref = gsl_matrix_get( LM -> alpha , i , j ) ;
      gsl_matrix_set( LM -> alpha_new , i , j , ref ) ;
      gsl_matrix_set( LM -> alpha_new , j , i , ref ) ;
    }
    gsl_matrix_set( LM -> alpha_new , i , i , 
		    ( 1.0 + Lambda ) * gsl_matrix_get( LM -> alpha , i , i ) ) ;
  }

  // solves alpha[p][q] * delta( a[q] ) = beta[p] for delta
#ifdef LMSVD
  if( gsl_linalg_SV_decomp( LM -> alpha_new , LM -> V ,
			    LM -> S , LM -> work ) != GSL_SUCCESS ) {
    fprintf( stderr , "[LM] SVD decomp failed\n" ) ;
    return !GSL_SUCCESS ;
  }
  if( gsl_linalg_SV_solve( LM -> alpha_new , LM -> V ,
			   LM -> S , LM -> beta , LM -> delta ) != GSL_SUCCESS ) {
    fprintf( stderr , "[LM] SVD decomp failed\n" ) ;
    return !GSL_SUCCESS ;
  }
#else
  int signum ;
  if( gsl_linalg_LU_decomp( LM -> alpha_new , LM -> perm , &signum ) != GSL_SUCCESS ) {
    fprintf( stderr , "[LM] LU decomp broken?\n" ) ;
    return !GSL_SUCCESS ;
  }
  if( gsl_linalg_LU_solve( LM -> alpha_new , LM -> perm , LM -> beta , LM -> delta ) != GSL_SUCCESS ) {
    fprintf( stderr , "[LM] LU solve broken?\n" ) ;
    return !GSL_SUCCESS ;
  }
#endif
  
  // update fitparams
  for( i = 0 ; i < f -> NPARAMS ; i++ ) {
    f -> fparams[i] = LM -> old_params[i] + gsl_vector_get( LM -> delta , i ) ;
    #ifdef VERBOSE
    printf( "[LM] NEW PARAMS :: %f \n" , f -> fparams[i] ) ;
    #endif
  }

  // compute the predicted chi^2?
  register double loc_sumd = 0.0 , loc_sumB = 0.0 ;
  for( i = 0 ; i < f -> NPARAMS ; i++ ) {
    register const double LMi = gsl_vector_get( LM -> delta , i ) ;
    loc_sumd += LMi*LMi ;
    register double sumj = 0.0 ; // matrix-vector product LM.alpha[i][j]*LM.delta[j] ;
    for( size_t j = 0 ; j < f -> NPARAMS ; j++ ) {
      sumj += gsl_matrix_get( LM -> alpha , i , j ) * gsl_vector_get( LM -> delta , j ) ;
    }
    loc_sumB +=  LMi*sumj ;
  }
  #ifdef VERBOSE
  fprintf( stdout , "[LM] F(x+d) = %e + %e \n" ,
	   loc_sumd , loc_sumB ) ;
  #endif
  
  LM -> pred = loc_sumd + loc_sumB/2. ;
  
  // compute new chisq
  fdesc.F( f -> f , data , f -> fparams ) ;
  return compute_chisq( *f , W , f -> CORRFIT ) ;
}

static void
init_LM( struct lmstep *LM ,
	 const struct ffunction f ,
	 const size_t Nlogic )
{
  LM->alpha     = gsl_matrix_alloc( Nlogic , Nlogic ) ;
  LM->alpha_new = gsl_matrix_alloc( Nlogic , Nlogic ) ;
  LM->beta      = gsl_vector_alloc( Nlogic ) ;
  LM->delta     = gsl_vector_alloc( Nlogic ) ;
  LM->perm      = gsl_permutation_alloc( Nlogic ) ;
  // allocate the y-data
  LM -> Nsum = f.N ;
  switch( f.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : LM -> Nsum = f.N * f.N ; break ;
  }
  LM -> y = malloc( LM -> Nsum * sizeof( double ) ) ;
#ifdef LMSVD
  LM->V = gsl_matrix_alloc( Nlogic , Nlogic ) ;
  LM->S = gsl_vector_alloc( Nlogic ) ;
  LM->work = gsl_vector_alloc( Nlogic ) ;
#endif
}

static void
free_LM( struct lmstep *LM )
{
  // free gsl stuff
  gsl_vector_free( LM -> beta ) ;
  gsl_matrix_free( LM -> alpha ) ;
  gsl_matrix_free( LM -> alpha_new ) ;
  gsl_permutation_free( LM -> perm ) ;
  gsl_vector_free( LM -> delta ) ;
  if( LM -> old_params != NULL ) {
    free( LM -> old_params ) ;
  }
  if( LM -> y != NULL ) {
    free( LM -> y ) ;
  }
#ifdef LMSVD
  gsl_matrix_free( LM -> V ) ;
  gsl_vector_free( LM -> work ) ;
  gsl_vector_free( LM -> S ) ;
#endif
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
  const size_t LMMAX = 5000 ; 

  double chisq_diff = 1E20 , Lambda = 1. ;
  size_t iters = 0 , i ;

  // lambda growth and shrinkage factors
  const double Dfac = 10 , Mfac = 4 ;

  // get priors
  Fit -> f.Prior = Fit -> Prior ;
  
  // allocate alpha, beta, delta and permutation matrices
  struct lmstep LM ;
  init_LM( &LM , Fit -> f , Fit -> Nlogic ) ;
  
  // set the old parameters
  //double old_params[ Fit -> Nlogic ] ;
  LM.old_params = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    LM.old_params[ i ] = Fit -> f.fparams[ i ] ;
  }

  // evaluate the function, its first and second derivatives
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
  #ifdef WITH_D2_DERIVS
  Fit -> d2F( Fit -> f.d2f , data , Fit -> f.fparams ) ;
  #endif
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

  // get alpha and beta for the new set of f
  get_alpha_beta( &LM , Fit -> f , W ) ;

  // loop until chisq evens out
  while( chisq_diff > TOL && iters < LMMAX ) {

    const double new_chisq = lm_step( &Fit -> f , &LM , *Fit , 
				      data , W , Lambda ) ;

    #ifdef VERBOSE
    fprintf( stdout , "[ML] chis :: %f %f %e \n" , new_chisq , 
	     Fit -> f.chisq , fabs( Fit -> f.chisq - new_chisq ) ) ;
    #endif
    
    // update lambda shrinking the size if we are close to a solution
    if( new_chisq <= Fit -> f.chisq ) {
      // update lambda
      Lambda /= Dfac ;
      chisq_diff = fabs( Fit -> f.chisq - new_chisq ) ;
      Fit -> f.chisq = new_chisq ;
      for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
	LM.old_params[ i ] = Fit -> f.fparams[ i ] ;
      }
      // update derivatives and alpha and beta
      Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
      #ifdef WITH_D2_DERIVS
      Fit -> d2F( Fit -> f.d2f , data , Fit -> f.fparams ) ;
      #endif
      get_alpha_beta( &LM , Fit -> f , W ) ;
    } else {
      for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
	Fit -> f.fparams[ i ] = LM.old_params[ i ] ;
      }
      Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ; // reset f
      Lambda *= Mfac ;
    }

    // if lambda misses a minimum tell us
    if( Lambda < 1E-32 || Lambda > 1E32 ) {
      fprintf( stderr , "[LM] Lambda is out of bounds %e \n" , Lambda ) ;
      iters = LMMAX ;
      break ;
    }
  
    // increment the number of iterations
    iters++ ;
  } // end of while loop

#ifdef VERBOSE
  // tell us how many iterations we hit
  if( iters == LMMAX ) {
    fprintf( stdout , "\n[LM] stopped by max iterations %zu -> Chidiff %e\n" ,
	     iters , chisq_diff ) ;
  } else {
    fprintf( stdout , "\n[LM] FINISHED in %zu iterations \n" , iters ) ;
  }

  fprintf( stdout , "[LM] chisq :: %e \n\n" , Fit -> f.chisq ) ;
  // tell us the fit parameters
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    fprintf( stdout , "PARAM_%zu :: %1.15e \n" , i , Fit -> f.fparams[i] ) ;
  }
#endif
  
  free_LM( &LM ) ;
  
  return iters ;
}
