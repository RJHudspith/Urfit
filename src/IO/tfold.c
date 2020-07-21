/**
   @file tfold.c
   @brief time correlation function folding
 */
#include <complex.h>

#include "gens.h"
#include "read_corr.h"

#define REAL

#ifdef REAL
static double part( const double complex C ) { return creal( C ) ; }
#elif (defined IMAG)
static double part( const double complex C ) { return cimag( C ) ; }
#elif (defined CABS)
static double part( const double complex C ) { return cabs( C ) ; }
#elif (defined RABS)
static double part( const double complex C ) { return fabs( creal( C ) ) ; }
#elif (defined IABS)
static double part( const double complex C ) { return fabs( cimag( C ) ) ; }
#else
static double part( const double complex C ) { return creal( C ) ; }
#endif

// map the input file selections to correlators in the file
double complex*
map_correlator( const struct traj Traj ,
		const char *str ,
		const double *mompoint ,
		const size_t Nlt )
{
  double complex *C = malloc( Nlt * sizeof( double complex ) ) ;
  size_t idx , Nsum = 0.0 ;
  
  // set correlator to zero
  for( idx = 0 ; idx < Nlt ; idx++ ) {
    C[idx] = 0.0 ;
  }
  
  const size_t *mapGs , *mapGk ;
  const size_t Vmap[3]   = { 0 , 1 , 2 } ;
  const size_t Amap[3]   = { 6 , 7 , 8 } ;
  const size_t Tijmap[3] = { 10 , 11 , 12 } ;
  const size_t Titmap[3] = { 13 , 14 , 15 } ;
  
  bool trigger1 = false , trigger2 = false ;

  if( Traj.Gs == Vi  ) { mapGs = Vmap ;   trigger1 = true ; }
  if( Traj.Gs == Ai  ) { mapGs = Amap ;   trigger1 = true ; }
  if( Traj.Gs == Tij ) { mapGs = Tijmap ; trigger1 = true ; }
  if( Traj.Gs == Tit ) { mapGs = Titmap ; trigger1 = true ; }
      
  if( Traj.Gk == Vi  ) { mapGk = Vmap ;   trigger2 = true ; }
  if( Traj.Gk == Ai  ) { mapGk = Amap ;   trigger2 = true ; }
  if( Traj.Gk == Tij ) { mapGk = Tijmap ; trigger2 = true ; }
  if( Traj.Gk == Tit ) { mapGk = Titmap ; trigger2 = true ; }
      
  if( trigger1 == true || trigger2 == true ) {
    for( idx = 0 ; idx < 3 ; idx++ ) {
      // TT
      if( trigger1 == true && trigger2 == true ) {
	if( get_correlator( C , str , mapGk[idx] , mapGs[idx] ,
			    mompoint , Nlt ) == FAILURE ) {
	  return NULL ;
	}
	Nsum++ ;
	// TF
      } else if( trigger1 == true && trigger2 == false ) {
	if( get_correlator( C , str , Traj.Gk , mapGs[idx] ,
			    mompoint , Nlt ) == FAILURE ) {
	  return NULL ;
	}
	Nsum++ ;
	// FT
      } else {
	if( get_correlator( C , str , mapGk[idx] , Traj.Gs ,
			    mompoint , Nlt ) == FAILURE ) {
	  return NULL ;
	}
	Nsum++ ;
      }
    }	  
  } else {
    if( get_correlator( C , str , Traj.Gk , Traj.Gs ,
			mompoint , Nlt ) == FAILURE ) {
      return NULL ;
    }
    Nsum++ ;
  }

  // normalize C by Nsum
  for( idx = 0 ; idx < Nlt ; idx++ ) {
    C[ idx ] /= Nsum ;
  }
  
  return C ;
}

// fold the correlator symmetrically
int
time_fold( struct resampled *sample ,
	   const double complex *C ,
	   const size_t LT ,
	   const fold_type fold ,
	   const size_t meas )
{
  const size_t L2 = LT/2+1 ;
  size_t t ;
  switch( fold ) {
  case PLUS_PLUS :
    sample[0].resampled[meas] = part( C[0] ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * part( C[t] + C[LT-t] ) ;
    }
    break ;
  case PLUS_MINUS :
    sample[0].resampled[meas] = part( C[0] ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = 0.5 * part( C[t] - C[LT-t] ) ;
    }
    break ;
  case MINUS_PLUS :
    sample[0].resampled[meas] = -part( C[0] ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = -0.5 * part( C[t] - C[LT-t] ) ;
    }
    break ;
  case MINUS_MINUS :
    sample[0].resampled[meas] = -part( C[0] ) ;
    for( t = 1 ; t < L2 ; t++ ) {
      sample[t].resampled[meas] = -0.5 * part( C[t] + C[LT-t] ) ;
    }
    break ;
  case NOFOLD :
    for( t = 0 ; t < LT ; t++ ) {
      sample[t].resampled[meas] = part( C[t] ) ;
    }
    break ;
  case NOFOLD_MINUS :
    for( t = 0 ; t < LT ; t++ ) {
      sample[t].resampled[meas] = -part( C[t] ) ;
    }
    break ;
  case TDER :
    for( t = 0 ; t < LT ; t++ ) {
      sample[t].resampled[meas] = 0.5*( part( C[(t+LT-1)%LT] + C[(t+1)%LT] ) ) ;
    }
    break ;
  case NOFOLD_SWAPT :
    sample[0].resampled[meas] = part( C[0] ) ;
    for( t = 1 ; t < LT ; t++ ) {
      sample[t].resampled[meas] = ( part( C[LT-t]) ) ;
    }
    break ;
  case NOFOLD_MINUS_SWAPT :
    sample[0].resampled[meas] = -part( C[0] ) ;
    for( t = 1 ; t < LT ; t++ ) {
      sample[t].resampled[meas] = -( part( C[LT-t]) ) ;
    }
    break ;
  }
  return SUCCESS ;
}
