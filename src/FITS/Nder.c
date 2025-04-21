/**
   @file Nder.c
   @brief numerical finite difference code
 */
#include "gens.h"

//#define VERBOSE

double
Nder( double (*f) ( const struct x_desc X ,
		    const double *fparams ,
		    const size_t Npars ) ,
      const struct x_desc X ,
      const size_t Npars ,
      const double *fparams ,
      const size_t idx ,        // parameter id
      const size_t Nparams )
{
  double Plus[ Nparams ] , Minus[ Nparams ] ;
  double h = 1E-2 , TOL = 1.0 , TOLPREV = UNINIT_FLAG , dfidx = 0 ;
  size_t inc = 0 ;
  
  // copy the temps
  memcpy( Plus  , fparams , Nparams * sizeof( double ) ) ;
  memcpy( Minus , fparams , Nparams * sizeof( double ) ) ;
  
  while( TOL > 1E-15 && inc < 200 ) {

    // set the plus and minus terms
    Plus[ idx ]  = fparams[ idx ] + h ;
    Minus[ idx ] = fparams[ idx ] - h ;

    #ifdef VERBOSE
    fprintf( stdout , "[NDER] PLUS %f MINUS %f \n"  ,
	     Plus[idx] , Minus[idx] ) ;
    #endif

    const double trial = ( f( X , Plus , Npars ) -
			   f( X , Minus , Npars ) ) / ( 2*h ) ;

    TOL = fabs( trial - dfidx ) ;

    #ifdef VERBOSE
    fprintf( stdout , "[NDER] TOL -> %e \n" , TOL ) ;
    #endif

    // if the solution is growing then break
    if( TOL > TOLPREV ) break ;

    h *= 0.1 ;
    dfidx = trial ;
    TOLPREV = TOL ;
    
    inc++ ;
  }

  #ifdef VERBOSE
  fprintf( stdout , "Finished :: dfidx %f \n" , dfidx ) ;
  #endif

  /*
  if( dfidx == UNINIT_FLAG || TOL > 1E-1 ) {
    fprintf( stderr , "[NDER] failed to converge :: %zu %e \n" , inc , TOL ) ;
    return sqrt(-1) ;
  }
  */
  
  return dfidx ;
}
    
