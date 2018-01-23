/**
   @file correlator.c
   @brief correlator analysis
 */
#include "gens.h"

#include "blackbox.h"
#include "decays.h"
#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"
#include "stats.h"
#include "write_flat.h"

//#define MATRIX_PRONY

void 
matrix_prony( const struct input_params *Input )
{
  const size_t Nstates = Input -> Fit.N ;
  size_t i , j , k , l , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    const size_t Ndata = Input -> Data.Ndata[i] ;
    double data[ Ndata ] ;
    double masses[ Nstates ][ Ndata ] ;

    struct resampled *effmass =
      malloc( Nstates * Ndata * sizeof( struct resampled ) ) ;
    for( j = 0 ; j < Nstates * Ndata ; j++ ) {
      effmass[j] = init_dist( NULL , Input -> Data.x[shift].NSAMPLES ,
			      Input -> Data.x[shift].restype ) ;
    }
    
    for( k = 0 ; k < Input -> Data.x[shift].NSAMPLES ; k++ ) {
      for( j = shift ; j < shift + Ndata ; j++ ) {
	data[ j - shift ] = Input -> Data.y[j].resampled[k] ;
      }
      blackbox( data , Ndata , Nstates , masses ) ;
      for( l = 0 ; l < Nstates ; l++ ) {
	for( j = 0 ; j < Ndata ; j++ ) {
	  effmass[j + l*Nstates].resampled[k] = masses[l][j] ;
	}
      }
      printf( "[PRONY] Samples %1.2f percent done\n" ,
	      100. * k/(double)Input->Data.x[shift].NSAMPLES ) ;
    }

    // same for the average
    for( j = shift ; j < shift + Ndata ; j++ ) {
      data[ j - shift ] = Input -> Data.y[j].avg ;
    }
    blackbox( data , Ndata , Nstates , masses ) ;
    for( l = 0 ; l < Nstates ; l++ ) {
      for( j = 0 ; j < Ndata ; j++ ) {
	effmass[j + l*Nstates].avg = masses[l][j] ;
      }
    } 
    
    for( j = 0 ; j < Nstates * Ndata ; j++ ) {
      compute_err( &effmass[j] ) ;
    }

    for( l = 0 ; l < Nstates ; l++ ) {
      for( j = 0 ; j < Ndata ; j++ ) {
	printf( "%zu %e %e \n" , j , effmass[j + l*Nstates].avg , effmass[j + l*Nstates].err ) ;
      }
      printf( "\n" ) ;
    }

    for( j = 0 ; j < Nstates * Ndata ; j++ ) {
      free( effmass[j].resampled ) ;
    }
    free( effmass ) ;
    
    shift += Ndata ;
  }

  return ;
}

int
correlator_analysis( struct input_params *Input )
{
  size_t i ;

#ifdef MATRIX_PRONY
  printf( "Matrix prony\n" ) ;

  matrix_prony( Input ) ;
#endif
  
  printf( "Effmass\n" ) ;
  
  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , LOG_EFFMASS ) ;

#ifdef FIT_EFFMASS
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    equate( &Input -> Data.y[i] , effmass[i] ) ;
  }
#endif

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;


  // compute the fractional error
  size_t shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t j ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      printf( "[FRAC] %e %e \n" , Input -> Data.x[j].avg ,
	      100 * ( Input -> Data.y[j].err / Input -> Data.y[j].avg ) ) ;
    }
    printf( "[FRAC] \n" ) ;
    shift += Input -> Data.Ndata[i] ;
  }

  // normalise to 1
  /*
  struct resampled div = init_dist( &Input -> Data.y[0] ,
				    Input -> Data.y[0].NSAMPLES ,
				    Input -> Data.y[0].restype ) ;
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t j ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      divide( &(Input -> Data.y[j]) , div ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  free( div.resampled ) ;
  */

  
  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  // compute a decay constant
  if( Input -> Fit.Fitdef == PP_AA_WW ||
      Input -> Fit.Fitdef == PP_AA ||
      Input -> Fit.Fitdef == PPAA ) {

    write_flat_dist( &Fit[0] , &Fit[0] , 1 , "Mass.flat" ) ;
    
    struct resampled dec = decay( Fit , *Input ) ;

    write_flat_dist( &dec , &dec , 1 , "Decay.flat" ) ;

    divide( &Fit[0] , dec ) ;

    raise( &Fit[0] , 2 ) ;
    
    printf( "(M/F)^2 %e %e \n" , Fit[0].avg , Fit[0].err ) ;

    write_flat_dist( &Fit[0] , &Fit[0] , 1 , "MovFsq.flat" ) ;
    
    free( dec.resampled ) ;
  }

  // write out a flat file
  if( Input -> Fit.Fitdef == EXP ) {
    //write_flat_single( &Fit[1] , "Mass.flat" ) ;
    struct resampled mpi2 = init_dist( NULL ,
				       Fit[1].NSAMPLES ,
				       Fit[1].restype ) ;
    equate_constant( &mpi2 , 0.088447949604 ,
		     Fit[1].NSAMPLES , Fit[1].restype ) ;
    write_flat_dist( &Fit[1] , &mpi2 , 1 , "Mass.flat" ) ;
    free( mpi2.resampled ) ;
  }
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
