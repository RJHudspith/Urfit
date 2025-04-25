/**
   minimize chi^2 by Powell's method
 */
#include "gens.h"

#include <limits.h>
#include <assert.h>

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h"

// minimises the parameters in p[n] copies over and returns chisq
int
powell_iter( void *fdesc ,
	     const void *data ,
	     const double **W ,
	     const double TOL )
{
  const int MAX_ITERS = 600 ;
  int iters = 0 ;

  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  const int n = Fit -> Nlogic ;
  
  // get priors and evaluate initial chisq and stuff
  Fit -> f.Prior = Fit -> Prior ;
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ; 
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
  
  // need this due to NR's crazy global variable solution
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , Fit -> f.N ) ;
  copy_ffunction( &f2 , Fit -> f ) ;
  
  double *fret = &Fit -> f.chisq ;
  
  // set the search directions to be the identity matrix
  double dirs[n][n] , p[n] , pt[n] ;
  for( int i = 0 ; i < n ; i++ ) {
    pt[i] = p[i] = Fit -> f.fparams[i] ;
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
      linmin( n , p , xit , fret , &f2 , fdesc , W , data ) ;
      if( (fptt - *fret) > del ) {
	del = fptt - *fret ;
	ibig = i ;
      }
    }
    // leave
    if( (2*(fp - *fret)) <= (TOL*(fabs(fp)+fabs(*fret))+1E-28 ) ) {  
      break ;
    }
    for( int j = 0 ; j < n ; j++ ) {
      ptt[j] = 2*p[j] - pt[j] ;
      xit[j] = p[j] - pt[j] ;
      pt[j]  = p[j] ;
      Fit -> f.fparams[j] = ptt[j] ;
    }
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
    fptt = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
    if( fptt < fp ) {
      const double t = 2*(fp-2*(*fret)+fptt)*sqrt(fp-*fret-del)-del*sqrt(fp-fptt) ;
      if( t < 0.0 ) {
	linmin( n , p , xit , fret , &f2 , fdesc , W , data ) ;
	for( int j = 0 ; j < n ; j++ ) {
	  dirs[ibig][j] = dirs[n-1][j] ;
	  dirs[n-1][j]  = xit[j] ;
	}
      }
    }
  }

  // gotta do this
  free_ffunction( &f2 , Fit -> Nlogic ) ;
    
  // here we reached the end and it is time to set everything
  for( int j = 0 ; j < n ; j++ ) {
    Fit -> f.fparams[j] = p[j] ;
  }
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = *fret ;
  
  return iters == MAX_ITERS ? FAILURE : iters ;
}
