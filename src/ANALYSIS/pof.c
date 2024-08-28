/**
   @file pof.c
   @brief compute the pencil of functions of a correlator (see hep-lat 1404.4029)
 */
#include "gens.h"

#include "effmass.h"
#include "gevp.h"
#include "stats.h"
#include "resampled_ops.h"
#include "fit_and_plot.h"

static const int spacing = 1;

// Jboxes 3,6

static const int t0 = 4 ; //2 ; //1*spacing ;
static const int tp = 8 ; //4 ; //2*spacing ;

static int
write_evalues( struct resampled *evalues ,
	       const size_t Ndata ,
	       const size_t N ,
	       const size_t t0 )
{
  size_t i , j = 0 , k ;
  for( i = 0 ; i < N ; i++ ) {
    char str[ 256 ] ;
    sprintf( str , "Evalue.%zu.flat" , i ) ;
    FILE *file = fopen( str , "w" ) ;

    fprintf( file , "%u\n" , evalues[ Ndata*i ].restype ) ;
    fprintf( file , "%zu\n" , Ndata ) ;
    
    for( j = 0 ; j < Ndata ; j++ ) {
      compute_err( &evalues[ j + Ndata*i ] ) ;
      printf( "%zu %e %e\n" , j ,
	      evalues[ j + Ndata*i ].avg ,
	      evalues[ j + Ndata*i ].err ) ;

      // write out the eigenvalues
      fprintf( file , "%zu\n" , evalues[ j + Ndata*i ].NSAMPLES ) ;
      for( k = 0 ; k < evalues[ j + Ndata*i ].NSAMPLES ; k++ ) {
	fprintf( file , "%f %1.15e\n" , (double)j-t0 , evalues[ j + Ndata*i ].resampled[k] ) ;
      }
      //
      fprintf( file , "AVG %f %1.15e\n" , (double)j-t0 , evalues[ j + Ndata*i ].avg ) ;
    }
    fclose( file ) ;
  }

  return SUCCESS ;
}

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
	const size_t tp = (t+spacing*n+spacing*m)%LT ;
	y[ idx ] = init_dist( &Input -> Data.y[tp] ,
			      Input -> Data.y[tp].NSAMPLES ,
			      Input -> Data.y[tp].restype ) ;
	printf( "y[%zu][%zu][%zu] = %e\n" , m , n , t , y[ idx ].avg ) ;
	idx++ ;
      }
    }
  }

  fprintf( stdout , "Solving GEVP\n" ) ;

  // compute evalues
  struct resampled *evalues = solve_GEVP( y ,
					  Input -> Data.Ndata[0] ,
					  Input -> Fit.N ,
					  Input -> Fit.M ,
					  t0 , tp ) ;

  // write them out
  write_evalues( evalues , Input -> Data.Ndata[0] , Input -> Fit.N , t0 ) ;
  
  // set input to the first evalue
  size_t i , j ;
  for( i = 0 ; i < 1 ; i++ ) {
    for( size_t j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
      equate( &Input -> Data.y[ j+i*Input->Data.Ndata[0] ] ,
	      evalues[ j + i*Input->Data.Ndata[0] ] ) ;
      printf( "Evalue_%zu %e %e \n" , i ,
	      Input -> Data.y[j+i*Input->Data.Ndata[0]].avg ,
	      Input -> Data.y[j+i*Input->Data.Ndata[0]].err ) ;
    }
    printf( "\n" ) ;
  }

  fprintf( stdout , "Effective mass\n") ;
  
  // compute an effective mass
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    //divide_constant( &Input -> Data.y[ i ] , 1E14 ) ;
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;
 
  for( i = 0 ; i < LT*(Input->Fit.M*Input->Fit.N) ; i++ ) {
    free( y[i].resampled ) ;
  }
  free( y ) ;

  fprintf( stdout , "Fit?\n") ;

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
