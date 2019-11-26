/**
   @file pof.c
   @brief compute the pencil of functions of a correlator (see hep-lat 1404.4029)
 */
#include "gens.h"

#include "effmass.h"
#include "gevp.h"
#include "resampled_ops.h"

// pencil of functions looks like
// C(t)   C(t+1)   .. C(t+N)
// C(t+2) C(t+2)   .. C(t+N+1)
// ....   ....     .. ....
// C(t+M) C(t+M+1) .. C(t+N+M)
//
// this looks a lot like the matrix prony method come to think about it
int
pof_analysis( struct input_params *Input )
{
  const int t0 = 1 ;
  struct resampled *y = malloc( Input -> Data.Ndata[0] * ( Input -> Fit.N*Input -> Fit.M )*sizeof( struct resampled ) ) ;

  const size_t LT = Input -> Data.Ndata[0] ;
  size_t t , m , n , idx = 0 ;
  for( m = 0 ; m < Input -> Fit.M ; m++ ) {
    for( n = 0 ; n < Input -> Fit.N ; n++ ) {
      for( t = 0 ; t < Input -> Data.Ndata[0] ; t++ ) {
	const size_t tp = (t+n+m)%LT ;
	y[ idx ] = init_dist( &Input -> Data.y[tp] ,
			      Input -> Data.y[tp].NSAMPLES ,
			      Input -> Data.y[tp].restype ) ;

	printf( "y[%zu][%zu][%zu] = %f\n" , m , n , t , y[ idx ].avg ) ;
	idx++ ;
      }
    }
  }

  // compute evalues
  struct resampled *evalues = solve_GEVP( y ,
					  Input -> Data.Ndata[0] ,
					  Input -> Fit.N ,
					  Input -> Fit.M ,
					  t0 , t0+1 ) ;

  // set input to the first evalue
  size_t i ;
  for( i = 0 ; i < Input -> Data.Ndata[0] ; i++ ) {
    equate( &Input -> Data.y[i] , evalues[i] ) ;
    printf( "TEST %e %e \n" , Input -> Data.y[i].avg ,
	    Input -> Data.y[i].err ) ;
  }

  // compute an effective mass
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;
 
  for( i = 0 ; i < LT*(Input->Fit.M*Input->Fit.N) ; i++ ) {
    free( y[i].resampled ) ;
  }
  free( y ) ;
  
  return SUCCESS ;
}
