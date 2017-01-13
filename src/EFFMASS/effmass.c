/**
   @file effmass.c
   @brief effective mass solutions
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "gens.h"
#include "fitfunc.h"
#include "resampled_ops.h"

#define EFFMASSTOL 1E-30

// set it to zero
static void
zero_effmass( struct resampled *effmass ,
	      const struct resampled bootavg )
{
  equate_constant( effmass , 0.0 , 
		   bootavg.NSAMPLES , 
		   bootavg.restype ) ;
  return ;
}

// acosh breaks down if acosh( n < 1 )
static bool
acosh_range( const struct resampled effmass )
{
  size_t i ;
  for( i = 0 ; i < effmass.NSAMPLES ; i++ ) {
    if( effmass.resampled[ i ] < 1 ) {
      return true ;
    }
  }
  return false ;
}

// make sure this is in the right range for the atanh
static bool
atanh_range( const struct resampled effmass ) 
{
  size_t i ;
  for( i = 0 ; i < effmass.NSAMPLES ; i++ ) {
    if( effmass.resampled[ i ] < -1 ||
	effmass.resampled[ i ] > 1 ) {
      return true ;
    }
  }
  return false ;
}

// if the data doesn't work we set the whole thing to zero
static bool
log_range( const struct resampled log )
{
  size_t i ;
  for( i = 0 ; i < log.NSAMPLES ; i++ ) {
    if( log.resampled[ i ] < 0 ) {
      return true ;
    }
  }
  return false ;
}

// there are some edge cases
static void
edge_cases( struct resampled *effmass ,
	    const struct resampled *bootavg ,
	    const size_t j ,
	    const size_t NDATA )
{
  if( j == 0 ) {
    equate( &effmass[j] , bootavg[j+1] ) ;
    divide( &effmass[j] , bootavg[j] ) ;
    if( log_range( effmass[j] ) ) {
      zero_effmass( &effmass[j] , bootavg[j] ) ;
    } else {
      res_log( &effmass[j] ) ;
      mult_constant( &effmass[j] , -1.0 ) ;
    }
  } else if( j == NDATA - 1 ) {
    equate( &effmass[j] , bootavg[j] ) ;
    divide( &effmass[j] , bootavg[j-1] ) ;
    // just in case there is a sign flip
    if( log_range( effmass[j] ) ) {
      return zero_effmass( &effmass[j] , bootavg[j] ) ;
    } else {
      res_log( &effmass[j] ) ;
      mult_constant( &effmass[j] , -1.0 ) ;
    }
  }
  return ;
}

// only works for evenly spaced separations atmo
struct resampled**
effective_mass( const struct resampled **bootavg ,
		const size_t *NDATA ,
		const size_t NSLICES ,
		const effmass_type type ,
		const size_t t0 )
{
  struct resampled **effmass = malloc( NSLICES * sizeof( struct resampled* ) ) ;
  size_t i ;
  #pragma omp parallel for private(i)
  for( i = 0 ; i < NSLICES ; i++ ) {
    effmass[ i ] = malloc( NDATA[i] * sizeof( struct resampled ) ) ;

    struct resampled temp ;
    size_t j ;

    switch( type ) {
    case LOG_EFFMASS :
      for( j = 0 ; j < NDATA[i] ; j++ ) {
	effmass[i][j].resampled = malloc( bootavg[i][j].NSAMPLES * sizeof( double ) ) ;
	if( log_range( bootavg[i][j] ) == true ) {
	  zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	} else {
	  // log effective mass?
	  if( j == 0 || j == NDATA[i] - 1 ) {
	    edge_cases( effmass[i] , bootavg[i] , j , NDATA[i] ) ;
	  } else {
	    equate( &effmass[i][j] , bootavg[i][j+1] ) ;
	    divide( &effmass[i][j] , bootavg[i][j] ) ;
	    // if it is still negative we set to zero
	    if( effmass[i][j].avg < EFFMASSTOL ) {
	      zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	    } else {
	      res_log( &effmass[i][j] ) ;
	    }
	    //
	    mult_constant( &effmass[i][j] , -1.0 ) ;
	  }
	}
	//
      }
      break ;
    case LOG2_EFFMASS :
      temp.resampled = malloc( bootavg[i][0].NSAMPLES * sizeof( double ) ) ;
      for( j = 0 ; j < NDATA[i] ; j++ ) {
	effmass[i][j].resampled = malloc( bootavg[i][j].NSAMPLES * sizeof( double ) ) ;
	if( log_range( bootavg[i][j] ) == true ) {
	  zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	} else {
	  // log effective mass?
	  if( j == 0 || j == NDATA[i] - 1 ) {
	    edge_cases( effmass[i] , bootavg[i] , j , NDATA[i] ) ;
	  } else {
	    equate( &effmass[i][j] , bootavg[i][j] ) ;
	    subtract( &effmass[i][j] , bootavg[i][j-1] ) ;
	    equate( &temp , bootavg[i][j+1] ) ;
	    subtract( &temp , bootavg[i][j] ) ;
	    divide( &effmass[i][j] , temp ) ;
	    // if it is still negative we set to zero
	    if( effmass[i][j].avg < EFFMASSTOL ) {
	      zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	    } else {
	      res_log( &effmass[i][j] ) ;
	    }
	    //
	  }
	}
      }
      free( temp.resampled ) ;
      break ;
      // cosh effective mass?
    case ACOSH_EFFMASS :
      for( j = 0 ; j < NDATA[i] ; j++ ) {
	effmass[i][j].resampled = malloc( bootavg[i][j].NSAMPLES * sizeof( double ) ) ;
	// set to zero if it is messed up
	if( acosh_range( bootavg[i][j] ) == true ) {
	  zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	} else {
	  if( j == 0 || j == NDATA[i] - 1 ) {
	    edge_cases( effmass[i] , bootavg[i] , j , NDATA[i] ) ;
	  } else {
	    equate( &effmass[i][j] , bootavg[i][j-1] ) ;
	    add( &effmass[i][j] , bootavg[i][j+1] ) ;
	    divide( &effmass[i][j] , bootavg[i][j] ) ;
	    mult_constant( &effmass[i][j] , 0.5 ) ;
	    // if it is in the wrong domain set to zero
	    if( acosh_range( effmass[i][j] ) ) {
	      zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	    } else {
	      res_acosh( &effmass[i][j] ) ;
	    }
	    //
	  }
	}
      }
      break ;
      // sinh effective mass
    case ASINH_EFFMASS :
      for( j = 0 ; j < NDATA[i] ; j++ ) {
	effmass[i][j].resampled = malloc( bootavg[i][j].NSAMPLES * sizeof( double ) ) ;
	if( fabs( bootavg[i][j].avg ) < EFFMASSTOL ) {
	  zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	} else {
	  if( j == 0 || j == NDATA[i] - 1 ) {
	    edge_cases( effmass[i] , bootavg[i] , j , NDATA[i] ) ;
	  } else {
	    equate( &effmass[i][j] , bootavg[i][j-1] ) ;
	    subtract( &effmass[i][j] , bootavg[i][j+1] ) ;
	    divide( &effmass[i][j] , bootavg[i][j] ) ;
	    mult_constant( &effmass[i][j] , 0.5 ) ;
	    // asinh works on all lengths
	    res_asinh( &effmass[i][j] ) ;
	  }
	}
      }
      break ;
    case ATANH_EFFMASS :
      temp.resampled = malloc( bootavg[i][0].NSAMPLES * sizeof( double ) ) ;
      // tanh effective mass?
      for( j = 0 ; j < NDATA[i] ; j++ ) {
	effmass[i][j].resampled = malloc( bootavg[i][j].NSAMPLES * sizeof( double ) ) ;
	if( atanh_range( bootavg[i][j] ) == true ) {
	  zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	} else {
	  if( j == 0 || j == NDATA[i] - 1 ) {
	    edge_cases( effmass[i] , bootavg[i] , j , NDATA[i] ) ;
	  } else {
	    // is E^{-m(t-1)} - E^{-m(t+1)}
	    equate( &effmass[i][j] , bootavg[i][j-1] ) ;
	    subtract( &effmass[i][j] , bootavg[i][j+1] ) ;
	    // is E^{-m(t-1)} + E^{-m(t+1)}
	    equate( &temp , bootavg[i][j-1] ) ;
	    add( &temp , bootavg[i][j+1] ) ;
	    // divide the two and take the tanh
	    divide( &effmass[i][j] , temp ) ;
	    // if it is still negative we set to zero
	    if( atanh_range( effmass[i][j] ) ) {
	      zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	    } else {
	      res_atanh( &effmass[i][j] ) ;
	    }
	    //
	  }
	}
      }
      free( temp.resampled ) ;
      break ;
      // eigenvalue effective mass
    case EVALUE_EFFMASS :
      for( j = 0 ; j < NDATA[i] ; j++ ) {
	effmass[i][j].resampled = malloc( bootavg[i][j].NSAMPLES * sizeof( double ) ) ;
	if( log_range( bootavg[i][j] ) == true ) {
	  zero_effmass( &effmass[i][j] , bootavg[i][j] ) ;
	} else {
	  equate( &effmass[i][j] , bootavg[i][j] ) ;
	  res_log( &effmass[i][j] ) ;
	  if( j > 0 ) {
	    mult_constant( &effmass[i][j] , -1.0/((double)j-t0) ) ;
	  }
	}
      }
      break ;
    }
  }
  return effmass ;
}

#undef EFFMASSTOL
