/**
   @file BFGS.c
   @brief BFGS minimizer for our \chi^2 functional
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h"

#include <assert.h>

//#define VERBOSE

#ifdef VERBOSE
// have a look at the hessian mayhaps
static void
printH( const int n , const double H[n][n] )
{
  printf( "\n H \n " ) ;
  for( int i = 0  ; i < n ; i++ ) {
    for( int j = 0 ; j < n ; j++ ) {
      printf( "%e " , H[i][j] ) ;
    }
    printf( "\n" ) ;
  }
}
#endif

// bfgs iterations
int
BFGS_iter( void *fdesc ,
	   const void *data ,
	   const double **W ,
	   const double TOL )
{
  // point to the fit descriptor
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  
  // iterations and CG maximum iterations
  size_t iters = 0 ;
  const size_t BFGSMAX = 8000 ;

  // allocate the fitfunction
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , Fit -> f.N ) ;
  copy_ffunction( &f2 , Fit -> f ) ;

  // get priors
  f2.Prior = Fit -> f.Prior = Fit -> Prior ;

  // inverse hessian estimate for initial guess set for now to be the identity matrix
  double H[ Fit -> Nlogic ][ Fit -> Nlogic ] ;  
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;

  // array temporaries of gradients and things
  double grad[ Fit -> Nlogic ] , s[ Fit -> Nlogic ] , p[ Fit -> Nlogic ] ;
  double y[ Fit -> Nlogic ] , Hy[ Fit -> Nlogic ] ;
  get_gradient( grad , W , Fit ) ;

  // set initial Hessian to be the identity
  for( int i = 0 ; i < Fit->Nlogic ; i++ ) {
    memset( H[i] , 0. , Fit->Nlogic*sizeof(double) ) ;
    H[i][i] = 1. ;
  }

  double chisq_diff = 10 , chiprev = 123456789 , chinew , alpha = 1 ;
  while( chisq_diff > TOL && iters < BFGSMAX ) {
    // update f which has the gradient direction in
    chinew = Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
    chisq_diff = fabs(chinew-chiprev) ;   
    chiprev = chinew ;
    // obtain direction P = -H_{ij}\Delta_j 
    for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
      register double sum = 0. ;
      for( int j = 0 ; j < Fit -> Nlogic ; j++ ) {
	sum += H[i][j]*grad[j] ;
      }
      p[i] = sum ;
    }
    // line search in this direction & set s = alpha*p and x = x + s
    alpha = line_search( &f2 , Fit -> f , p , *Fit , data , W ) ;
    for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
      s[i] = alpha*p[i] ;
      Fit -> f.fparams[i] += s[i] ;
    }
    // recompute gradients and put in y = \Del
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
    Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
    // get new descent direction (-grad hence all the fucking signs)
    get_gradient( y , W , Fit ) ;

    // update estimate for inverse Hessian
    double yHy = 0.0 , sty = 0 , sumy = 0 , sums = 0. ;
    for( int i = 0  ; i < Fit -> Nlogic ; i++ ) {
      double tmp = (y[i] - grad[i]) ; grad[i] = y[i] ; y[i] = tmp ;
      sty  += -s[i]*y[i] ;
      sumy += sqrt( y[i]*y[i] ) ;
      sums += sqrt( s[i]*s[i] ) ;
      register double sum = 0 ;
      for( int j = 0 ; j < Fit -> Nlogic ; j++ ) {
	sum += H[i][j]*y[j] ;
      }
      Hy[i] = -sum ;
      yHy += y[i]*sum ;
    }

    // these are normalisers so if they are super small we are fucked and we just leave
    if( fabs(yHy) < 1E-15 || fabs(sty) < 1E-15 ) {
      break ;
    }
    
    // NR values
    const double fad = 1./yHy , fac = 1./sty , fae = yHy ;
    if( fac > sqrt(1E-15*sumy*sums) ) {
      // this is the magic vector that gets added to make it a full BFGS
      for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
	y[i] = fac*s[i] - fad*Hy[i] ;
      }
      // cool outerproducts s.sT and such
      for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
	const register double A =  fac*s[i] ;
	const register double B = -fad*Hy[i] ;
	const register double C =  fae*y[i] ;
	// try and urge the compiler to FMA these
	for( int j = 0 ; j < Fit -> Nlogic ; j++ ) {
	  H[i][j] += A*s[j] ;
	  H[i][j] += B*Hy[j] ;
	  H[i][j] += C*y[j] ;
	}
      }
    }    
    iters++ ;
  }
  // tell us how many iterations we hit
#ifdef VERBOSE
  if( iters == CGMAX ) {
    printf( "\n[BFGS] stopped by max iterations %zu \n" , iters ) ;
  } else {
    printf( "\n[BFGS] FINISHED in %zu iterations \n" , iters ) ;
  }
  printf( "[BFGS] chisq :: %e | Diff -> %e \n\n" , Fit -> f.chisq , chisq_diff ) ;
  for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
    printf( "PARAMS :: %e \n" , Fit -> f.fparams[i] ) ;
  }
#endif
  free_ffunction( &f2 , Fit -> Nlogic ) ; 
  return iters ;
}

// clean up the defines
#ifdef VERBOSE
  #undef VERBOSE
#endif

