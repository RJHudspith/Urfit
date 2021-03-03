/**
   @file pof.c
   @brief compute the pencil of functions of a correlator (see hep-lat 1404.4029)
 */
#include "gens.h"

#include "effmass.h"
#include "gevp.h"
#include "resampled_ops.h"
#include "fit_and_plot.h"

// pencil of functions looks like
// C(t)   C(t+1)   .. C(t+N)
// C(t+1) C(t+2)   .. C(t+N+1)
// ....   ....     .. ....
// C(t+M) C(t+M+1) .. C(t+N+M)
//
// this looks a lot like the matrix prony method come to think about it
int
pof_analysis( struct input_params *Input )
{
  const size_t N = Input -> Fit.N ;
  
  const int t0 = 1 ;
  struct resampled *y = malloc( Input -> Data.Ndata[0] * ( Input -> Fit.N*Input -> Fit.M )*sizeof( struct resampled ) ) ;

  size_t t , m , n , idx = 0 ;
  /*
  for( t = 0 ; t < Input -> Data.Ndata[0] ; t++ ) {
    divide_constant( &Input -> Data.y[t] , 1E16 ) ;
  }
  */

  const size_t LT = Input -> Data.Ndata[0] ;
  for( m = 0 ; m < Input -> Fit.M ; m++ ) {
    for( n = 0 ; n < Input -> Fit.N ; n++ ) {
      for( t = 0 ; t < Input -> Data.Ndata[0] ; t++ ) {
	const size_t tp = (t+n+m)%LT ;
	y[ idx ] = init_dist( &Input -> Data.y[tp] ,
			      Input -> Data.y[tp].NSAMPLES ,
			      Input -> Data.y[tp].restype ) ;

	printf( "y[%zu][%zu][%zu] = %e\n" , m , n , t , y[ idx ].avg ) ;
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
  size_t i , j ;
  for( i = 0 ; i < 1; i++ ) {
    for( size_t j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
      equate( &Input -> Data.y[ j+i*Input->Data.Ndata[0] ] ,
	      evalues[ j + i*Input->Data.Ndata[0] ] ) ;
      printf( "TEST %e %e \n" ,
	      Input -> Data.y[j+i*Input->Data.Ndata[0]].avg ,
	      Input -> Data.y[j+i*Input->Data.Ndata[0]].err ) ;
    }
  }

  // compute an effective mass
  struct resampled *effmass = effective_mass( Input , LOG_EFFMASS ) ;

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    divide_constant( &Input -> Data.y[ i ] , 1E14 ) ;
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;
 
  for( i = 0 ; i < LT*(Input->Fit.M*Input->Fit.N) ; i++ ) {
    free( y[i].resampled ) ;
  }
  free( y ) ;

  // perform a fit ?
  const size_t Nsim_prev = Input -> Data.Nsim , Ntot_prev = Input -> Data.Ntot ;
  Input -> Data.Nsim = N ;
  Input -> Data.Ntot = N*Input->Data.Ndata[0] ;
  Input -> Fit.N = 1 ;
  Input -> Fit.M = 1 ;

  fprintf( stdout , "Dofit\n" ) ;
  
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  Input -> Data.Nsim = Nsim_prev ;
  Input -> Data.Ntot = Ntot_prev ;
  
  return SUCCESS ;
}


// pencil of functions looks like
// C(t)   C(t+1)   .. C(t+N)
// C(t+2) C(t+2)   .. C(t+N+1)
// ....   ....     .. ....
// C(t+M) C(t+M+1) .. C(t+N+M)
//
// this looks a lot like the matrix prony method come to think about it
int
pof_analysis_fixed( struct input_params *Input )
{
  const size_t N = Input -> Fit.N ;
  
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

  struct resampled *evalues = solve_GEVP_fixed( y ,
						Input -> Data.Ndata[0] ,
						Input -> Fit.N ,
						Input -> Fit.M ,
						t0 , t0+1 ) ;

  size_t i , j ;
  for( size_t j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    equate( &Input -> Data.y[ j ] , evalues[ j ] ) ;
    res_log( &Input -> Data.y[ j ] ) ;
    divide_constant( &Input -> Data.y[ j ] , t0 ) ;
    printf( "TEST %e %e \n" ,
	    Input -> Data.y[j].avg ,
	    Input -> Data.y[j].err ) ;
  }
 
  for( i = 0 ; i < LT*(Input->Fit.M*Input->Fit.N) ; i++ ) {
    free( y[i].resampled ) ;
  }
  free( y ) ;

  // perform a fit ?
  const size_t Nsim_prev = Input -> Data.Nsim , Ntot_prev = Input -> Data.Ntot ;
  Input -> Data.Nsim = N ;
  Input -> Data.Ntot = N*Input->Data.Ndata[0] ;
  Input -> Fit.N = 0 ;
  Input -> Fit.M = 1 ;

  fprintf( stdout , "Dofit\n" ) ;
  
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  Input -> Data.Nsim = Nsim_prev ;
  Input -> Data.Ntot = Ntot_prev ;
  
  return SUCCESS ;
}
