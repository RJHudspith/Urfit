/**
   conjugate gradient routine
 */
#include "gens.h"

#include "chisq.h"
#include "ffunction.h"
#include "line_search.h"
#include "summation.h"

//#define VERBOSE

// cg iteration, this is pretty good now
int
cg_iter( void *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // point to the fit descriptor
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
  
  // iterations and CG maximum iterations
  size_t iters = 0 ;
  const size_t CGMAX = 8000 ;

  // allocate the fitfunction
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , Fit -> f.N ) ;
  copy_ffunction( &f2 , Fit -> f ) ;

  // get priors
  f2.Prior = Fit -> f.Prior = Fit -> Prior ;

  // evaluate the function, and its first derivatives
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

  // allocate conjugate directions
  double *s      = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  double *old_df = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  get_gradient( old_df , W , Fit ) ;
  memcpy( s , old_df , Fit->Nlogic*sizeof(double) ) ;

  // line search the SD step
  double alpha = line_search( &f2 , Fit -> f , old_df , s , *Fit , data , W ) ;
  for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
    Fit -> f.fparams[i] += alpha*s[i] ;
  }
  
  double chisq_diff = 10 , chiprev = 123456789 , chinew ;
  while( sqrt(chisq_diff) > TOL && iters < CGMAX ) {
    // update f which has the gradient direction in
    Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
    Fit -> dF( Fit -> f.df , data , Fit -> f.fparams ) ;
    chinew = Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;
    chisq_diff = fabs( chinew - chiprev ) ;
    chiprev = chinew ;
    // compute beta using polyak - ribiere
    register double num = 0.0 , denom = 0.0 ;
    double newdf[ Fit -> Nlogic ] ;
    get_gradient( newdf , W , Fit ) ;
    for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
      num   += newdf[i]*(newdf[i]-old_df[i]) ;
      //num   += newdf[i]*newdf[i] ; Fletcher reeves
      denom += old_df[i]*old_df[i] ;
      old_df[i] = newdf[i] ;
    }
    double beta = fmax( 0. , num/denom ) ;
    for( int i = 0 ; i < Fit->Nlogic ; i++ ) {
      s[i] = newdf[i] + beta*( s[i] ) ;
    }
    // update with a line search
    alpha = line_search( &f2 , Fit -> f , old_df , s , *Fit , data , W ) ;
    for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
      Fit -> f.fparams[i] += alpha*s[i] ;
    }
    iters++ ;
  }
  // tell us how many iterations we hit
#ifdef VERBOSE
  if( iters == CGMAX ) {
    printf( "\n[CG] stopped by max iterations %zu \n" , iters ) ;
  } else {
    printf( "\n[CG] FINISHED in %zu iterations \n" , iters ) ;
  }
  printf( "[CG] chisq :: %e | Diff -> %e \n\n" , Fit -> f.chisq , chisq_diff ) ;
  for( int i = 0 ; i < Fit -> Nlogic ; i++ ) {
    printf( "PARAMS :: %e \n" , Fit -> f.fparams[i] ) ;
  }
#endif
  // free the directions
  if( s != NULL ) {
    free( s ) ;
  }
  if( old_df != NULL ) {
    free( old_df ) ;
  }
  return iters ;
}

// clean up the defines
#ifdef VERBOSE
  #undef VERBOSE
#endif

