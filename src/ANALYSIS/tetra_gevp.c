/**
   @file correlator.c
   @brief correlator analysis
 */
#include "gens.h"

#include "decays.h"
#include "effmass.h"
#include "fit_and_plot.h"
#include "gevp.h"
#include "init.h"
#include "resampled_ops.h"
#include "stats.h"

//#define FIT_EFFMASS
//#define COMBINE

int
tetra_gevp_analysis( struct input_params *Input )
{
  const size_t t0 = 3 ;
  const size_t N = (size_t)sqrt( Input -> Data.Nsim ) ;
  const size_t LT = Input -> Traj[0].Dimensions[3]/2 ;
  size_t i , j , k ;

  printf( "[TET GEVP] In here :: %zu %zu %zu\n" , t0 , N , LT ) ;
  
  // make a correlator matrix
  // C_ij / sqrt( C_ii * C_jj )
  size_t shift = 0 ;
  fprintf( stdout , "[TET GEVP] correlator matrix\n" ) ;
  for( i = 0 ; i < N ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      printf( "%e " , Input -> Data.y[ t0 + LT*( j +  N*i ) ].avg /
	      sqrt( Input -> Data.y[ t0 + LT*( j + N*j ) ].avg *
		    Input -> Data.y[ t0 + LT*( i + N*i ) ].avg ) ) ;
    }
    printf( "\n" ) ;
  }

  // subtract the shift automatically
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    subtract_constant( &Input -> Data.x[ i ] , (double)t0 ) ;
  }
  
  // allocate the eigenvalues
  struct resampled *evalues = malloc( Input -> Data.Ndata[0] * N *
				      sizeof( struct resampled ) ) ;
  for( j = 0 ; j < Input ->Data.Ndata[0]*N ; j++ ) {
    evalues[j].resampled = malloc( Input -> Data.y[0].NSAMPLES *
				   sizeof( double ) ) ;
    evalues[j].restype = Input -> Data.y[0].restype ;
    evalues[j].NSAMPLES = Input -> Data.y[0].NSAMPLES ;
  }

  // ugh loop order is all weird
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {

    // temporary matrices C0 and C1 and vector re_evalues
    double C0[ Input -> Data.Nsim ] ;
    double C1[ Input -> Data.Nsim ] ;
    double re_evalues[ N ] ;
    
    for( k = 0 ; k < Input -> Data.y[0].NSAMPLES ; k++ ) {

      size_t shift = 0 ;
      for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
	C0[i] = Input -> Data.y[ j  + shift ].resampled[k] ;
	C1[i] = Input -> Data.y[ t0 + shift ].resampled[k] ;
	shift += Input -> Data.Ndata[i] ;
      }

      // compute eigenvalues
      if( solve_gevp( re_evalues , C0 , C1 , N , j<t0 , false ) == FAILURE ) {
        // maybe do something?
      }

      // poke into solution
      for( i = 0 ; i < N ; i++ ) {
	evalues[j+Input->Data.Ndata[0]*i].resampled[ k ] = re_evalues[ i ] ;
      }
      //
    }

    // redo for the average
    size_t shift = 0 ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      C0[i] = Input -> Data.y[ j  + shift ].avg ;
      C1[i] = Input -> Data.y[ t0 + shift ].avg ;
      shift += Input -> Data.Ndata[i] ;
    }

    // compute eigenvalues
    if( solve_gevp( re_evalues , C0 , C1 , N , j<t0 , true ) == FAILURE ) {
      // maybe do something?
    }

    // poke into solution
    for( i = 0 ; i < N ; i++ ) {
      evalues[j+Input->Data.Ndata[0]*i].avg = re_evalues[ i ] ;
    }
    //
  }

  for( i = 0 ; i < N ; i++ ) {
    char str[ 256 ] ;
    sprintf( str , "Evalue.%zu.flat" , i ) ;
    FILE *file = fopen( str , "w" ) ;

    fprintf( file , "%zu\n" , evalues[ j + Input->Data.Ndata[0]*i ].restype ) ;
    fprintf( file , "%zu\n" , Input->Data.Ndata[0] ) ;
    
    for( j = 0 ; j < Input->Data.Ndata[0] ; j++ ) {
      compute_err( &evalues[ j + Input->Data.Ndata[0]*i ] ) ;
      printf( "%zu %e %e\n" , j ,
	      evalues[ j + Input->Data.Ndata[0]*i ].avg ,
	      evalues[ j + Input->Data.Ndata[0]*i ].err ) ;

      // write out the eigenvalues
      fprintf( file , "%zu\n" , evalues[ j + Input->Data.Ndata[0]*i ].NSAMPLES ) ;
      for( k = 0 ; k < evalues[ j + Input->Data.Ndata[0]*i ].NSAMPLES ; k++ ) {
	fprintf( file , "%f %1.15e\n" , (double)j , evalues[ j + Input->Data.Ndata[0]*i ].resampled[k] ) ;
      }
      //
      fprintf( file , "AVG %f %1.15e\n" , (double)j , evalues[ j + Input->Data.Ndata[0]*i ].avg ) ;
    }
    fclose( file ) ;
  }


  // copy over to Input
  for( i = 0 ; i < N ; i++ ) {
    for( j = 0 ; j < Input->Data.Ndata[0] ; j++ ) {
      equate( &Input -> Data.y[j+i*Input->Data.Ndata[0]] ,
	      evalues[ j + i*Input->Data.Ndata[0] ] ) ;
    }
    printf( "\n" ) ;
  }

  // free the eigenvalues
  for( i = 0 ; i < Input ->Data.Ndata[0] * N  ; i++ ) {
    free( evalues[i].resampled ) ;
  }
  free( evalues ) ;

  const size_t Nsim_prev = Input -> Data.Nsim , Ntot_prev = Input -> Data.Ntot ;
  Input -> Data.Nsim = N ;
  Input -> Data.Ntot = N*Input->Data.Ndata[0] ;
  
  // compute an effective mass
  struct resampled *effmass = effective_mass( Input , ACOSH_EFFMASS ) ;

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  Input -> Data.Nsim = 1 ;
  Input -> Data.Ntot = Input->Data.Ndata[0] ;
  Input -> Fit.Nlogic = 2 ;

  // fit the ground state
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  Input -> Data.Nsim = Nsim_prev ;
  Input -> Data.Ntot = Ntot_prev ;

  return SUCCESS ;
}
