/**
   minimize chi^2 by Powell's method
 */
#include "gens.h"

#include <limits.h>
#include <assert.h>

#include "chisq.h"
#include "ffunction.h"

#define MAXNPAR (24)

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
  loc_fdesc.F( loc_f.f , loc_data , loc_f.fparams ) ;
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
  while( fabs( x3 - x0 ) > tol*(fabs(x1)+fabs(x2) ) ) {
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
      fprintf( stderr , "Precision %e\n" , fabs(x3-x0) ) ;
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

// need to convert here between our fit function and simple 1D function for the line search
static void
linmin( const int n ,
	double p[n] ,
	double xi[n] ,
	double *fret ,
	const struct fit_descriptor *fdesc ,
	const double **W ,
	const void *data )  
{
  // probably not going to have more parameters than this
  assert( n < 24 ) ;
  copy_ffunction( &loc_f , fdesc -> f ) ;
  loc_fdesc = *fdesc ;
  ncom = fdesc -> Nlogic ; // # of logical parameters
  loc_W = W ;
  loc_data = data ;
  memcpy( pcom  , p  , n*sizeof( double ) ) ;
  memcpy( xicom , xi , n*sizeof( double ) ) ;
  double a = 0 , b = 5 , c = 10 , xmin = 0 ;
  bracket( f1dim , &a , &b , &c ) ;  
  *fret = line_NR( f1dim , a , b , c , &xmin , 1E-3 ) ;  
  for( int j = 0 ; j < n ; j++ ) {
    xi[j] *= xmin ;
    p[j] += xi[j] ;
  }
}

// minimises the parameters in p[n] copies over and returns chisq
int
powell_iter( void *fdesc ,
	     const void *data ,
	     const double **W ,
	     const double TOL )
{
  const int MAX_ITERS = 200 ;
  int iters = 0 ;
  // need this due to NR's crazy global variable solution
#pragma omp critical
  {  
    struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
    const int n = Fit -> Nlogic ;
  
    struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , Fit -> f.N ) ;  
    loc_f = f2 ; 
  
    // get priors and evaluate initial chisq and stuff
    Fit -> f.Prior = Fit -> Prior ;
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ; 
    Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
  
    double *fret = &Fit -> f.chisq ;
  
    // set the search directions to be the identity matrix
    double dirs[n][n] , p[n] , pt[n] ;
    for( int i = 0 ; i < n ; i++ ) {
      pt[i] = p[ i ] = Fit -> f.fparams[i] ;
      for( int j = 0 ; j < n ; j++ ) {
	dirs[i][j] = (i==j)? 1.0 : 0.0 ;
      }
    }
    double xit[n] , ptt[n] , fptt = 0 , fp = 0 ;
  
    // value of the function at point "p"
    for( iters = 0 ; iters < MAX_ITERS ; iters++ ) {
      double del = 0. ;
      int ibig = 0 ;
      fp = *fret ;
      for( int i = 0 ; i < n ; i++ ) {
	memcpy( xit , dirs[i] , n*sizeof(double)) ;
	fptt = *fret ;
	linmin( n , p , xit , fret , fdesc , W , data ) ;
	if( (fptt - *fret) > del ) {
	  del = fptt - *fret ;
	  ibig = i ;
	}
      }
      // leave
      if( (2*(fp - *fret)) <= (TOL*(fabs(fp)+fabs(*fret))+1E-28 ) ) {  
	break ;
      }
      //
      for( int j = 0 ; j < n ; j++ ) {
	ptt[j] = 2*p[j] - pt[j] ;
	xit[j] = p[j] - pt[j] ;
	pt[j]  = p[j] ;
      }
      for( int j = 0 ; j < n ; j++ ) {
	Fit -> f.fparams[j] = ptt[j] ;
      }
      Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
      fptt = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
      if( fptt < fp ) {
	const double t = 2*(fp-2*(*fret)+fptt)*sqrt(fp-*fret-del)-del*sqrt(fp-fptt) ;
	if( t < 0.0 ) {
	  linmin( n , p , xit , fret , fdesc , W , data ) ;
	  for( int j = 0 ; j < n ; j++ ) {
	    dirs[ibig][j] = dirs[n-1][j] ;
	    dirs[n-1][j]  = xit[j] ;
	  }
	}
      }
    }
    // free the fitfunction
    free_ffunction( &f2 , Fit -> Nlogic ) ;
  
    // here we reached the end and it is time to set everything
    for( int j = 0 ; j < n ; j++ ) {
      Fit -> f.fparams[j] = p[j] ;
    }
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
    Fit -> f.chisq = *fret ;
  }
  return iters == MAX_ITERS ? FAILURE : iters ;
}
