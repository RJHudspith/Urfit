/**
   @file correlator.c
   @brief correlator analysis
 */
#include "gens.h"

#include "decays.h"
#include "effmass.h"
#include "fit_and_plot.h"
#include "gevp.h"
#include "make_xmgrace.h"
#include "init.h"
#include "resampled_ops.h"
#include "stats.h"

static const int t0 = 2 ;


//#define FIT_EFFMASS
//#define COMBINE

static int
write_fitmass_graph( FILE *file , 
		     const struct resampled mass ,
		     const double lo ,
		     const double hi ,
		     const int t0 )
{
  fprintf( file , "%e %e\n%e %e\n\n" , lo+t0 , mass.err_hi , hi+t0 , mass.err_hi ) ; 
  fprintf( file , "%e %e\n%e %e\n\n" , lo+t0 , mass.avg , hi+t0 , mass.avg ) ;
  fprintf( file , "%e %e\n%e %e\n\n" , lo+t0 , mass.err_lo , hi+t0 , mass.err_lo ) ;
  return SUCCESS ;
}

static int
write_evalues( struct resampled *evalues ,
	       const size_t Ndata ,
	       const size_t N ,
	       const size_t t0 )
{
  size_t i , j , k ;
  for( i = 0 ; i < N ; i++ ) {
    char str[ 256 ] ;
    sprintf( str , "Evalue.%zu.flat" , i ) ;
    FILE *file = fopen( str , "w" ) ;

    fprintf( file , "%u\n" , evalues[ j + Ndata*i ].restype ) ;
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

int
tetra_gevp_analysis( struct input_params *Input )
{
  //const size_t t0 = 2 ;
  const size_t N = Input -> Fit.N ;
  const size_t LT = Input -> Data.Ndata[0] ;
  size_t i , j , shift = 0 ; 

  printf( "[TET GEVP] In here :: %zu %zu %zu\n" , t0 , N , LT ) ;
  
  // make a correlator matrix
  // C_ij / sqrt( C_ii * C_jj )
  fprintf( stdout , "[TET GEVP] correlator matrix\n" ) ;
  size_t t ;
  for( t = 0 ; t < LT ; t++ ) {
    printf( "\n C(%zu) \n" , t ) ;
    for( i = 0 ; i < Input -> Fit.M ; i++ ) {
      for( j = 0 ; j < N ; j++ ) {
	printf( "%g " , Input -> Data.y[ t + LT*( j +  N*i ) ].avg
		/sqrt( Input -> Data.y[ t + LT*( j + N*j ) ].avg * Input -> Data.y[ t + LT*( i + N*i ) ].avg )
		) ;
      }
      printf( "\n" ) ;
    }
  }
  
  printf( "Divide ?\n" ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    struct resampled V = init_dist( &Input -> Data.y[i*Input->Data.Ndata[0]] ,
				    Input -> Data.y[i*Input->Data.Ndata[0]].NSAMPLES ,
				    Input -> Data.y[i*Input->Data.Ndata[0]].restype ) ;
				    
    for( j = 0 ; j < Input->Data.Ndata[0] ; j++ ) {
      //divide( &Input -> Data.y[j+i*Input->Data.Ndata[0]] , V ) ;
    }
    free( V.resampled ) ;
  }

  if( ( Input -> Fit.N*Input-> Fit.M ) != Input -> Data.Nsim ) {
    fprintf( stderr , "GEVP N,M not equal to Nsim\n" ) ;
    goto end ;
  }
  
  // compute evalues
  struct resampled *evalues = solve_GEVP( Input -> Data.y ,
					  Input -> Data.Ndata[0] ,
					  Input -> Fit.N ,
					  Input -> Fit.M ,
					  t0 ) ;
  
  if( evalues == NULL ) {
    goto end ;
  }

  // square root it if we need to
  if( Input -> Fit.N != Input -> Fit.M ) {
    for( i = 0 ; i < Input -> Fit.N ; i++ ) {
      for( j = 0 ; j < Input->Data.Ndata[0] ; j++ ) {
	root( &evalues[ j + i*Input->Data.Ndata[0] ] ) ;
      }
    }
  }
  
  // write them out
  write_evalues( evalues , Input -> Data.Ndata[0] , Input -> Fit.N , t0 ) ;

  for( i = 0 ; i < Input -> Fit.N ; i++ ) {
    for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
      equate( &Input -> Data.y[j+i*Input->Data.Ndata[0]] , evalues[ j + i*Input->Data.Ndata[0] ] ) ;
    }
  }

  // free the eigenvalues
  for( i = 0 ; i < Input ->Data.Ndata[0] * Input -> Fit.N  ; i++ ) {
    free( evalues[i].resampled ) ;
  }
  free( evalues ) ;

  const size_t Nsim_prev = Input -> Data.Nsim , Ntot_prev = Input -> Data.Ntot ;
  Input -> Data.Nsim = N ;
  Input -> Data.Ntot = N*Input->Data.Ndata[0] ;
  
  // compute an effective mass
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

  // subtract the ground state?
  for( j = 0 ; j < Input->Data.Ndata[0] ; j++ ) {
    struct resampled sub = init_dist( &effmass[j] ,
				      effmass[j].NSAMPLES ,
				      effmass[j].restype ) ;
    printf( "SUB %zu\n" , j ) ;
    for( i = 0 ; i < Input -> Fit.N ; i++ ) {
      subtract( &effmass[j+i*Input->Data.Ndata[0]] , sub ) ;
    }
    free( sub.resampled ) ;
  }

  for( i = 1 ; i < Input -> Fit.N ; i++ ) {
    for( j = 0 ; j < Input->Data.Ndata[0] ; j++ ) {
      printf( "Subbed %e %e %e \n" , Input -> Data.x[j].avg ,
	      effmass[j+i*Input->Data.Ndata[i]].avg ,
	      effmass[j+i*Input->Data.Ndata[i]].err ) ;
    }
    printf( "\n" ) ;
  }

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // subtract the shift automatically for the fit?
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    subtract_constant( &Input -> Data.x[ i ] , (double)t0 ) ;
  }

  Input -> Data.Nsim = 1 ;
  Input -> Data.Ntot = Input->Data.Ndata[0] ;
  Input -> Fit.N = Input -> Fit.M = 1 ;
  Input -> Fit.Nlogic = 2 ;

  
  // fit the ground state
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  // write out the fit result for reasons
  if( Input -> Fit.Fitdef == EXP ||
      Input -> Fit.Fitdef == COSH ) {
    FILE *massfile = fopen( "massfits.dat" , "w" ) ;
    size_t j ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      for( j = 0 ; j < 2*Input -> Fit.N ; j+= 2 ) {
	write_fitmass_graph( massfile , Fit[j+1] ,
			     Input -> Traj[i].Fit_Low ,
			     Input -> Traj[i].Fit_High ,
			     t0 ) ;
      }
      shift += Input -> Data.Ndata[i] ;
    }
    fclose( massfile ) ;
  }
  
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  Input -> Data.Nsim = Nsim_prev ;
  Input -> Data.Ntot = Ntot_prev ;

 end :

  return SUCCESS ;
}
