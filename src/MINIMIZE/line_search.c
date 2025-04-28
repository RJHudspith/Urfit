/**
   @file line_search.c

   various line search codes from numerical recipes (changed a lot because their code is incomprehensible and didn't parallelise) and myself
 */
#include "gens.h"

#include <limits.h>
#include <assert.h>
#include <omp.h>

#include "chisq.h"
#include "ffunction.h"
#include "summation.h"

// yeah Brent's method is a little faster
#define BRENT

// store all this crap
typedef struct {
  const int n ;
  double *p , *xi ;
  struct ffunction *f2 ;
  const struct fit_descriptor *fdesc ;
  const double **W ;
  const void *data ;
} ltemps ;

// hmmm this gets called a lot maybe I should make it faster
static inline double
ntrialf( ltemps tmp , const double x )
{
  for( int j = 0 ; j < tmp.n ; j++ ) {
    tmp.f2 -> fparams[j] = tmp.p[j] + x*tmp.xi[j] ;
  }
  tmp.fdesc -> F( tmp.f2 -> f , tmp.data , tmp.f2 -> fparams ) ;
  return compute_chisq( *tmp.f2 , tmp.W , tmp.f2 -> CORRFIT ) ; 
}

// compute the bracketing using the golden ratio from NR
static void
nbracket( ltemps tmp , double *a , double *b , double *c )
{ 
  const double GOLD = 1.618034 ;
  double ax = *a , bx = *b ;
  double fa = ntrialf( tmp , ax ) , fb = ntrialf( tmp , bx ) ;
  // swaps
  if( fb > fa ) {
    double tmp = ax ; ax = bx ; bx = tmp ;
    tmp = fa ; fa = fb ; fb = tmp ;
  }
  // guess for c
  double cx = bx + GOLD*(bx-ax) ;
  double fc = ntrialf( tmp , cx ) ;
  while( fb > fc ) {
    const double r = (bx-ax)*(fb-fc) , q = (bx-cx)*(fb-fa) ;
    double u = bx-((bx-cx)*q-(bx-ax)*r)/
      (2.*copysign(fmax(fabs(q-r),1E-20),q-r)) ;
    double ulim = bx + 111*(cx-bx) , fu = 1 ;
    if( (bx-u)*(u-cx) > 0 ) {
      fu = ntrialf( tmp , u ) ; 
      if( fu < fc ) {
	ax = bx ; bx = u  ;
	fa = fb ; fb = fu ;
	break ;
      } else if( fu > fb ) {
	cx = u  ; fc = fu ;
	break ;
      }
      u = cx + GOLD*(cx-bx) ;
      fu = ntrialf( tmp , u ) ; 
    } else if( (cx-u)*(u-ulim) > 0.0 ) {
      fu = ntrialf( tmp , u ) ; 
      if( fu < fc ) {
	bx = cx ; cx = u  ; u = u+GOLD*(u-cx) ;
	fb = fc ; fc = fu ; fu = ntrialf( tmp , u ) ;
      }
    } else if( (u-ulim)*(ulim-cx) >= 0.0 ) {
      u = ulim ;
      fu = ntrialf( tmp , u ) ; 
    } else {
      u = cx + GOLD*(cx-bx) ;
      fu = ntrialf( tmp , u ) ;
    }
    ax = bx ; bx = cx ; cx = u ;
    fa = fb ; fb = fc ; fc = fu ;
  }
  // logically a and c are our intervals and b lies between s.t. f(b) is lower
  *a = ax ; *b = bx ; *c = cx ;
}

#ifdef BRENT
// returns the function evaluated at xmin and xmin itself using Brent's method
static double
nline_NR( ltemps tmp , const double ax , const double bx , const double cx , double *xmin , const double tol )
{
  const int MAXLINE = 200 ;
  const double R = 0.6180339887498949 , C = 1.-R ;
  //const double ax = brac[0] , bx = brac[1] , cx = brac[2] ;
  double a = ax < cx ? ax : cx , b = ax > cx ? ax : cx ;
  double x = bx , w = bx , v = bx , e = 0.0 ;
  double fw , fv , fx , tol1 ;
  fw = fv = fx = ntrialf( tmp , x ) ;
  for( int iters = 0 ; iters < MAXLINE ; iters++ ) {
    const double xm = 0.5*(a+b) ;
    const double tol2 = 2*(tol1=tol*fabs(x)+1E-20) ;
    double d = 0. , u = 0 ;
    if( fabs(x-xm) <= (tol2-0.5*(b-a)) ) {
      break ;
    }
    // logic block for the parabolic step 
    if( fabs(e) > tol1 ) {
      const double r = (x-w)*(fx-fv) ;
      double q = (x-v)*(fx-fw) ;
      double p = (x-v)*q-(x-w)*r ;
      q = 2*(q-r) ;
      if( q > 0. ) p = -p ;
      q = fabs(q) ;
      double etemp = e ;
      // can the code ever get here on first iter as d is unitialized in NR
      e = d ;
      if( fabs(p) >= fabs( 0.5*q*etemp ) ||
	  p <= q*(a-x) ||
	  p >= q*(b-x) ) {
	d = C*(e=(x>xm?a-x:b-x)) ;
      } else {
	d = p/q ;
	u = x+d ;
	if( (u-a)<tol2 || (b-u)<tol2 ) {
	  d = copysign(tol1,xm-x) ;
	}
      }
    } else {
      d = C*(e=(x>xm?a-x:b-x)) ;
    }
    // single function evaluation
    u = (fabs(d) >= tol1 ? x+d : x+copysign(tol1,d)) ;
    const double fu = ntrialf( tmp , u  ) ;
    if( fu <= fx ) {
      if( u >= x ) a = x ; else b = x ;
      v = w ; w = x ; x = u ;
      fv = fw ; fw = fx ; fx = fu ;
    } else {
      if( u < x ) a = u ; else b = u ;
      if( fu <= fw || w == x ) {
	v  =  w ; w  =  u ;
	fv = fw ; fw = fu ;
      } else if( fu <= fv || v == x || v == w ) {
	v  =  u ;
	fv = fu ;
      }
    }
  }
  *xmin = x ;
  return fx ;
}
#else
// returns the function evaluated at xmin and xmin itself using golden section search
static double
nline_NR( ltemps tmp , const double a , const double b , const double c , double *xmin , const double tol )
{
  const int MAXLINE = 200 ;
  const double R = 0.6180339887498949 , C = 1.-R ;
  double x0 = a , x1 , x2 , x3 = c ;
  if( fabs( c - b ) > fabs( b - a ) ) {
    x1 = b ; x2 = b + C*(c-b) ; 
  } else {
    x2 = b ; x1 = b - C*(b-a) ;
  }
  double fx1 = ntrialf( tmp , x1 ) , fx2 = ntrialf( tmp , x2 ) ; 
  int iters = 0 ;
  while( fabs( x3 - x0 ) > tol*(fabs(x1)+fabs(x2))+1E-25 ) {
    if( fx2 < fx1 ) {
      x0 = x1 ; x1 = x2 ; x2 = R*x1 + C*x3 ; 
      fx1 = fx2 ; fx2 = ntrialf( tmp , x2 ) ; 
    } else {
      x3 = x2 ; x2 = x1 ; x1 = R*x2 + C*x0 ;
      fx2 = fx1 ; fx1 = ntrialf( tmp , x1 ) ; 
    }
    if( iters > MAXLINE ) {
      fprintf( stderr , "Line search failed to meet desired precision in %d iterations\n" , iters ) ;
      fprintf( stderr , "Precision %e\n" , fabs(x3-x0) ) ;
      break ;
    }
    iters++ ;
  }
  if( fx1 < fx2 ) {
    *xmin = x1 ; return fx1 ;
  } else {
    *xmin = x2 ; return fx2 ;
  }
}
#endif

// function is local and can be called in a thread-parallel fashion for Powell
double
linmin( const int n , double p[n] , double xi[n] ,
	double *fret , struct ffunction *f2 ,
	const struct fit_descriptor *fdesc ,
	const double **W ,
	const void *data )  
{
  ltemps tmp = { .n = n , .p = p , .xi = xi , .f2 = f2 ,
		 .fdesc = fdesc , .W = W , .data = data } ;
  double a = 0 , b = 0.5 , c = 1, xmin = 0 ;  
  nbracket( tmp , &a , &b , &c ) ;  
  *fret = nline_NR( tmp , a , b , c , &xmin , 1E-3 ) ;
  for( int j = 0 ; j < n ; j++ ) {
    xi[j] *= xmin ;
    p[j] += xi[j] ;
  }
  return xmin ;
}

// perform an step in the descent direction
static double
test_step( struct ffunction *f2 , 
	   const double *grad ,
	   const struct ffunction f1 ,
	   const struct fit_descriptor fdesc ,
	   const void *data ,
	   const double **W ,
	   const double alpha )
{
  // copy f1 to f2
  copy_ffunction( f2 , f1 ) ;
  for( int i = 0 ; i < fdesc.Nlogic ; i++ ) {
    f2 -> fparams[ i ] += alpha * grad[i] ;
  }
  fdesc.F( f2 -> f , data , f2 -> fparams ) ;
  return compute_chisq( *f2 , W , f2 -> CORRFIT ) ;
}

// backtracking line search and then golden ratio improvement
double
line_search( struct ffunction *f2 ,
	     const struct ffunction f1 ,
	     const double *grad ,
	     const double *descent , 
	     const struct fit_descriptor fdesc ,
	     const void *data ,
	     const double **W ,
	     double atrial )
{
  const double fac = 0.1 ;
  // perform a "backtracking line search" could use brent's method
  // do a rough search 100 -> 1E-16 for abest step 10
  double min = 123456789 ;
  double abest = 1E-15 ;

  const double f0 = test_step( f2 , descent , f1 , fdesc , data , W , 0.0 ) ;
  double armijo = 0 ;
  for( int i = 0 ; i < fdesc.Nlogic ; i++ ) {
    armijo += descent[i]*descent[i] ;
  }
  armijo *= fac ;
  
  while( atrial > 1E-16 ) {
    atrial *= fac ;
    double trial = test_step( f2 , descent , f1 , fdesc , data , W , atrial ) ;
    if( isnan( trial ) ) continue ;
    if( isinf( trial ) ) continue ;
    if( trial < min ) {
      abest = atrial ;
      min = trial ;
    }
    // backtrack until we first satisfy the armijo constraint
    if( trial < f0 - atrial*armijo ) break ;
  }

  // and then Brent/golden ratio search to refine the step
  const double R = 0.6180339887498949 ;
#ifdef BRENT
  const double C = 1.-R , tol = 1E-6 ;
  const int maxiters = 100 ;
  double ax = abest*fac , bx = abest , cx = abest/fac ;
  double a = ax < cx ? ax : cx , b = ax > cx ? ax : cx ;
  double x = bx , w = bx , v = bx , e = 0.0 ;
  double fw , fv , fx , tol1 ;
  fw = fv = fx = test_step( f2 , descent , f1 , fdesc , data , W , x ) ;
  int iters ;
  for( iters = 0 ; iters < maxiters ; iters++ ) {
    const double xm = 0.5*(a+b) ;
    const double tol2 = 2*(tol1=tol*fabs(x)+1E-20) ;
    double d = 0. , u = 0 ;
    if( fabs(x-xm) <= (tol2-0.5*(b-a)) ) {
      break ;
    }
    // logic block for the parabolic step 
    if( fabs(e) > tol1 ) {
      const double r = (x-w)*(fx-fv) ;
      double q = (x-v)*(fx-fw) ;
      double p = (x-v)*q-(x-w)*r ;
      q = 2*(q-r) ;
      if( q > 0. ) p = -p ;
      q = fabs(q) ;
      double etemp = e ;
      // can the code ever get here on first iter as d is unitialized in NR
      e = d ;
      if( fabs(p) >= fabs( 0.5*q*etemp ) ||
	  p <= q*(a-x) ||
	  p >= q*(b-x) ) {
	d = C*(e=(x>xm?a-x:b-x)) ;
      } else {
	d = p/q ;
	u = x+d ;
	if( (u-a)<tol2 || (b-u)<tol2 ) {
	  d = copysign(tol1,xm-x) ;
	}
      }
    } else {
      d = C*(e=(x>xm?a-x:b-x)) ;
    }
    // single function evaluation
    u = (fabs(d) >= tol1 ? x+d : x+copysign(tol1,d)) ;
    const double fu = test_step( f2 , descent , f1 , fdesc , data , W , u ) ;
    if( fu <= fx ) {
      if( u >= x ) a = x ; else b = x ;
      v = w ; w = x ; x = u ;
      fv = fw ; fw = fx ; fx = fu ;
    } else {
      if( u < x ) a = u ; else b = u ;
      if( fu <= fw || w == x ) {
	v  =  w ; w  =  u ;
	fv = fw ; fw = fu ;
      } else if( fu <= fv || v == x || v == w ) {
	v  =  u ;
	fv = fu ;
      }
    }
  }
  //printf( "Iters %d \n" , iters ) ;
  return x ;
#else
  double a = abest/fac , b = abest*fac ; //abest*fac , b = abest/fac ;
  double c = b - ( b-a )*R;
  double d = a + ( b-a )*R ;
  double fc = test_step( f2 , descent , f1 , fdesc , data , W , c ) ;
  double fd = test_step( f2 , descent , f1 , fdesc , data , W , d ) ;
  while( fabs( b-a ) > 1E-7*(fabs(b)+fabs(a)) ) {
    if( fc < fd ) {
      b = d ; d = c ; c = b - R*(b-a) ;
      fd = fc ; fc = test_step( f2 , descent , f1 , fdesc , data , W , c ) ;
    } else {
      a = c ; c = d ; d = a + R*(b-a) ;
      fc = fd ; fd = test_step( f2 , descent , f1 , fdesc , data , W , d ) ;
    }
  }
  return 0.5*(a+b) ;
#endif
}

// gets minus the derivative of the \chi^2 function i.e. the descent direction
void
get_gradient( double *grad ,
	      const double **W ,
	      const struct fit_descriptor *Fit )
{
  size_t Nsum = Fit -> f.N ;
  switch( Fit -> f.CORRFIT ) {
  case UNWEIGHTED : case UNCORRELATED : break ;
  case CORRELATED : Nsum = Fit -> f.N * Fit -> f.N ; break ;
  }
  double y[ Nsum ] ; 
  for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
    switch( Fit -> f.CORRFIT ) {
    case UNWEIGHTED :
      for( int j = 0 ; j < Fit -> f.N ; j++ ) {
	y[j] = -Fit -> f.df[i][j] * Fit -> f.f[j] ; 
      }
      break ;
    case UNCORRELATED :
      for( int j = 0 ; j < Fit -> f.N ; j++ ) {
	y[j] = -Fit -> f.df[i][j] * W[0][j] * Fit -> f.f[j] ; 
      }
      break ;
    case CORRELATED :
      for( int j = 0 ; j < Fit -> f.N ; j++ ) {
	for( int k = 0 ; k < Fit -> f.N ; k++ ) {
	  y[k+Fit->f.N*j] = -Fit -> f.df[i][j] * W[j][k] * Fit -> f.f[k] ; 
	}
      }
      break ;
    }
    grad[i] = kahan_summation( y , Nsum ) ;
    if( Fit -> f.Prior[i].Initialised == true ) {
      grad[i] -= ( Fit -> f.fparams[i] - Fit -> f.Prior[i].Val ) / 
	( Fit -> f.Prior[i].Err * Fit -> f.Prior[i].Err ) ;
    }
  }
}

#ifdef BRENT
  #undef BRENT
#endif
