/**
   @file simplex.c
   @bried minimization using the Nelder-Mead downhill simplex

   Also, fuck Numerical Recipes as I couldn't get their algorithm to make sense
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>

typedef struct {
  double *p ;
  double feval ;
} simplex ;

int compare( const void *a , const void *b )
{
  const simplex *A = (const simplex*)a , *B = (const simplex*)b ;
  return A -> feval > B -> feval ;
}

// average over n vertices omitting n+1 very importantly
static inline void
ave( const int n , const simplex s[n+1]  , double ave[n] ) {
  for( int j = 0 ; j < n ; j++ ) {
    register double sum = 0. ;
    for( int i = 0 ; i < n ; i++ ) {
      sum += s[i].p[j] ;
    }
    ave[j] = sum/n ;
  }
}

static double
evaluate_chisq( const int n ,
		const double x[n] ,
		struct ffunction *f2 ,
		const struct fit_descriptor *Fit ,
		const double **W ,
		const double *data )
{
  for( int i = 0 ; i < n ; i++ ) {
    f2 -> fparams[i] = x[i] ;
  }
  Fit -> F( f2 -> f , data , f2 -> fparams ) ;
  return compute_chisq( *f2 , W , f2 -> CORRFIT ) ;
}

// replace a simplex point
static inline void
replace( const int n , simplex *s , const double x[n] , const double fxn )
{
  memcpy( s->p , x , n*sizeof(double) ) ; s -> feval = fxn ;
}

// returns function eval at trial x
// xupdate is x = x_1 + multiplier*( x2 - x3 )
static inline double
trialx( const int n , double xnew[n] , const double multiplier ,
	const double x1[n] , const double x2[n] , const double x3[n] ,
	struct ffunction *f2 ,
	const struct fit_descriptor *Fit ,
	const double **W ,
	const double *data )
{
  for( int j = 0 ; j < n ; j++ ) {
    xnew[j] = x1[j] + multiplier*( x2[j] - x3[j] ) ;
  }
  return evaluate_chisq( n , xnew , f2 , Fit , W , data ) ;
}

int simplex_iter( void *fdesc ,
		  const void *data ,
		  const double **W ,
		  const double TOL )
{
  // default parameters
  const double alpha = 1 ;   // reflection
  const double beta = 2 ;    // expansion
  const double gamma = 0.5 ; // contraction
  const double delta = 0.5 ; // shrink
  const int MAX_ITERS = 5000 ;
  
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  const int n = Fit -> Nlogic ;
  
  // get priors and evaluate initial chisq and stuff
  Fit -> f.Prior = Fit -> Prior ;
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ; 
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
  
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , Fit -> f.N ) ;
  copy_ffunction( &f2 , Fit -> f ) ;
  
  simplex s[ n+1 ] ;
  s[0].p = malloc( n*sizeof(double) ) ;
  memcpy( s[0].p , Fit->f.fparams , n*sizeof(double) ) ;
  s[0].feval = evaluate_chisq( n , s[0].p , &f2 , Fit , W , data ) ;
  
  // probably there are smarter choices for these in the wild
  for( int i = 1 ; i < n+1 ; i++ ) {
    s[i].p = malloc( n*sizeof(double) ) ;
    for( int j = 0 ; j < n ; j++ ) {
      if( (i-1) == j ) {
	s[i].p[j] = (1.0+0.025)*s[0].p[j] ;
      } else {
	s[i].p[j] = 0.005 ;
      }
    }
    s[i].feval = evaluate_chisq( n , s[i].p , &f2 , Fit , W , data ) ;
  }

  int iters ;
  for( iters = 0 ; iters < MAX_ITERS ; iters++ ) {
    // sort to make sure best simplex is [0] and worst is [n]
    qsort( s , n+1 , sizeof(simplex) , compare ) ; 
    // convergence criteria
    double conv_v = 0 , conv_f = 0 ;
    for( int i = 1 ; i < n+1 ; i++ ) {
      double vec = 0 ;
      for( int j = 0 ; j < n ; j++ ) {
	vec += pow( s[i].p[j] - s[0].p[j] , 2 ) ;
      }
      conv_v = fmax( conv_v , sqrt(vec/n) ) ;
      conv_f = fmax( conv_f , sqrt(fabs( s[i].feval - s[0].feval )) ) ;
    }
    if( conv_v < TOL ) {
      break ;
    }
    // reflection step
    double xbar[ n ] , xr[ n ] , xo[ n ] , xe[ n ] , xi[ n ] ;
    ave( n , s , xbar ) ;
    const double fr = trialx( n , xr , alpha , xbar , xbar , s[n].p , &f2 , Fit , W , data ) ;    
    bool doshrink = false ;
    // expansion
    if( fr < s[0].feval ) {
      double xe[ n ] ;
      const double fe = trialx( n , xe , beta , xbar , xr , xbar , &f2 , Fit , W , data ) ;
      if( fe < fr ) {
	replace( n , &s[n] , xe , fe ) ;
      } else {
	replace( n , &s[n] , xr , fr ) ;
      }
    } else { // s[n-1].feval > fxr >= s[0].feval 
      if( fr < s[n-1].feval ) {
	replace( n , &s[n] , xr , fr ) ;
      } else { // fxr 
	if( fr < s[n].feval ) {
	  // outside contraction
	  double xo[ n ] ;
	  const double fo = trialx( n , xo , gamma , xbar , xr , xbar , &f2 , Fit , W , data ) ;
	  if( fo <= fr ) {
	    replace( n , &s[n] , xo , fo ) ;
	  } else {
	    doshrink = true ;
	  }
	} else {
	  // inside contraction
	  double xi[ n ] ;
	  const double fi = trialx( n , xi , -gamma , xbar , xr , xbar , &f2 , Fit , W , data ) ;
	  if( fi < s[n].feval ) {
	    replace( n , &s[n] , xi , fi ) ;
	  } else {
	    doshrink = true ;
	  }
	}
      }
    }
    // shrink if we have to
    if( doshrink == true ) {
      for( int i = 1 ; i < n+1 ; i++ ) {
	s[i].feval = trialx( n , s[i].p , delta , s[0].p , s[i].p , s[0].p , &f2 , Fit , W , data ) ;
      }
    }
  }
  
  // copy back to guess vector x and set function return pointer fret
  memcpy( Fit -> f.fparams , s[0].p , n*sizeof(double) ) ;
  Fit -> f.chisq = s[0].feval ;
  
  // free S[i].p here and we are done
  for( int i = 0 ; i < n+1 ; i++ ) {
    free( s[i].p ) ;
  }
  // gotta do this
  free_ffunction( &f2 , Fit -> Nlogic ) ;

  return iters == MAX_ITERS ? FAILURE : iters ;
}
