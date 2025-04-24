/**
   @file line_search.c
   @brief numerical recipes line search stuff
 */
#include "gens.h"

#include <limits.h>
#include <assert.h>

#include "chisq.h"
#include "ffunction.h"
#include "summation.h"

// put all the line search stuff in here
#define MAXNPAR (24)

//#define VERBOSE

// don't see much improvement between Golden Search and Brent's method honestly
//#define BRENT

// eww NR use some gross global functions this will be VERY problematic in parallel
static int ncom = 0 ;

// change with the number of threads...
static struct ffunction loc_f ;
static struct fit_descriptor loc_fdesc ;

// don't change
static const double **loc_W ;
static const void *loc_data ;

static double pcom[MAXNPAR] , xicom[MAXNPAR] ;

// gross conversion of our chisq to a one-dimensional function for the minimizers
static inline double
f1dim( const double x )
{
  // adjust the fit parameters
  for( int j = 0 ; j < ncom ; j++ ) {
    loc_f.fparams[j] = pcom[j] + x*xicom[j] ;
  }
  // compute the residual
  loc_fdesc.F(  loc_f.f  , loc_data , loc_f.fparams ) ;
  //loc_fdesc.dF( loc_f.df , loc_data , loc_f.fparams ) ;
  // return the chi^2
  return compute_chisq( loc_f , loc_W , loc_f.CORRFIT ) ;
}

// one dimensional function bracket using the golden ratio from NR
static void
bracket( double (*f)( const double x ) ,
	 double *a , double *b , double *c )
{
  const double GOLD = 1.618034 ;
  double ax = *a , bx = *b ;
  double fa = f(ax) , fb = f(bx) ;
  // swaps
  if( fb > fa ) {
    double tmp = ax ; ax = bx ; bx = tmp ;
    tmp = fa ; fa = fb ; fb = tmp ;
  }
  // guess for c
  double cx = bx + GOLD*(bx-ax) ;
  double fc = f(cx) , fu = 1 ;  
  while( fb > fc ) {
    const double r = (bx-ax)*(fb-fc) ;
    const double q = (bx-cx)*(fb-fa) ;
    double u = bx-((bx-cx)*q-(bx-ax)*r)/
      (2.*copysign(fmax(fabs(q-r),1E-20),q-r)) ;
    double ulim = bx + 111*(cx-bx) ;
    if( (bx-u)*(u-cx) > 0 ) {
      fu = f(u) ;
      if( fu < fc ) {
	ax = bx ;
	bx = u ;
	fa = fb ;
	fb = fu ;
	break ;
      } else if( fu > fb ) {
	cx = u ;
	fc = fu ;
	break ;
      }
      u = cx + GOLD*(cx-bx) ;
      fu = f(u) ;
    } else if( (cx-u)*(u-ulim) > 0.0 ) {
      fu = f(u) ;
      if( fu < fc ) {
	bx = cx ; cx = u  ; u = u+GOLD*(u-cx) ;
	fb = fc ; fc = fu ; fu = f(u) ;
      }
    } else if( (u-ulim)*(ulim-cx) >= 0.0 ) {
      u = ulim ;
      fu = f(u) ;
    } else {
      u = cx + GOLD*(cx-bx) ;
      fu = f(u) ;
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
line_NR( double (*f)( const double x ) ,
	  const double ax , const double bx , const double cx ,
	  double *xmin ,
	  const double tol )
{
  const int MAXLINE = 200 ;
  const double R = 0.6180339887498949 , C = 1.-R ;

  double a = ax < cx ? ax : cx ;
  double b = ax > cx ? ax : cx ;
  double x = bx , w = bx , v = bx , e = 0.0 ;
  double fw , fv , fx , tol1 ;
  fw = fv = fx = f( x ) ;
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
      // can the code ever get here as d is unitialized in NR
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
    const double fu = f(u) ;
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
// returns the function evaluated at xmin and xmin itself assumes a < b < c
static double
line_NR( double (*f)( const double x ) ,
	 const double a , const double b , const double c ,
	 double *xmin ,
	 const double tol )
{
  const int MAXLINE = 200 ;
  const double R = 0.6180339887498949 , C = 1.-R ;
  double x0 = a , x1 , x2 , x3 = c ;
  if( fabs( c - b ) > fabs( b - a ) ) {
    x1 = b ; x2 = b + C*(c-b) ; 
  } else {
    x2 = b ; x1 = b - C*(b-a) ;
  }
  double f1 = f(x1) , f2 = f(x2) ;
  int iters = 0 ;
  while( fabs( x3 - x0 ) > tol*(fabs(x1)+fabs(x2)+1E-25 ) ) {
    if( f2 < f1 ) {
      x0 = x1 ; x1 = x2 ; x2 = R*x1 + C*x3 ; 
      f1 = f2 ; f2 = f(x2) ;
    } else {
      x3 = x2 ; x2 = x1 ; x1 = R*x2 + C*x0 ;
      f2 = f1 ; f1 = f(x1) ;
    }
    #ifdef VERBOSE
    printf( "Line %e %e\n" , fabs(x3-x0) , tol*(fabs(x1)+fabs(x2))) ;
    #endif
    if( iters > MAXLINE ) {
      fprintf( stderr , "Line search failed to meet desired precision in %d iterations\n" , iters ) ;
      fprintf( stderr , "Precision %e | target %e\n" , fabs(x3-x0) , tol*(fabs(x1)+fabs(x2)) ) ;
      break ;
    }
    iters++ ;
  }
  #ifdef VERBOSE
  printf( "Line search convergence in %d iterations\n" , iters ) ;
  #endif
  if( f1 < f2 ) {
    *xmin = x1 ;
    return f1 ;
  } else {
    *xmin = x2 ;
    return f2 ;
  }
}
#endif

// need to convert here between our fit function and simple 1D function for the line search
void
linmin( const int n ,
	double p[n] ,
	double xi[n] ,
	double *fret ,
	struct ffunction *f2 , 
	const struct fit_descriptor *fdesc ,
	const double **W ,
	const void *data )  
{
  loc_f = *f2 ;
  // probably not going to have more parameters than this
  assert( n < 24 ) ;
  copy_ffunction( &loc_f , fdesc -> f ) ;
  loc_fdesc = *fdesc ;
  ncom = fdesc -> Nlogic ; // # of logical parameters
  loc_W = W ;
  loc_data = data ;
  memcpy( pcom  , p  , n*sizeof( double ) ) ;
  memcpy( xicom , xi , n*sizeof( double ) ) ;
  double a = 0 , b = 2 , c = 5 , xmin = 0 ;
  bracket( f1dim , &a , &b , &c ) ;  
  *fret = line_NR( f1dim , a , b , c , &xmin , 1E-3 ) ;
  for( int j = 0 ; j < n ; j++ ) {
    xi[j] *= xmin ;
    p[j] += xi[j] ;
  }
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
	     const double **W )
{
  const double fac = 0.1 ;

  // perform a "backtracking line search" should use brent's method
  // really and probably will do at some point, uses a golden ratio
  // search to find approximately the minimum
  // do a rough search 1^12 -> 1E-12 for abest step 100
  double min = 123456789 ;
  double atrial = 100 , abest = 1E-15 ;
  size_t iters = 0 ;
  while( atrial > 1E-16 ) {
    atrial *= fac ;
    double trial = test_step( f2 , descent , f1 , fdesc , data , W , atrial ) ;
    iters++ ;
    if( isnan( trial ) ) continue ;
    if( isinf( trial ) ) continue ;
    if( trial < min
	// Armijo
	// && trial > fx - atrial*0.1*descent[jidx]*descent[jidx]
	) {
      abest = atrial ;
      min = trial ;
    }
  }  
  
  double amid = abest ;
  double a = amid / fac , b = amid * fac ;
  // perform golden ratio line search
  const double R = 0.6180339887498949 ;
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
  #ifdef VERBOSE
  if( iters == 1000 ) {
    fprintf( stderr , "[LINE SEARCH] backtracking failed %e \n" , abest ) ;
  }
  printf( "[LINE SEARCH] abest :: %e \n" , abest ) ;
  #endif
  return (a+b)/2 ;
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
