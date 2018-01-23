/**
   @file Qcorr.c
   @brief topological correlator measurement
 */
#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "momenta.h"
#include "resampled_ops.h"

//#define SUBZERO

// the the topological corr < Q >^2 / V
int
fit_Qcorr( struct input_params *Input )
{   
  size_t i , j , shift = 0 ;
  const double fac = 1.0/ ( 10*10*10*10. ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift = j ;
  }
  
  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}

// partially sum the correlator \sum_{x=0}^{r^2} < q(x+y) q(x) >
int
fit_Qsusc( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // subtract the zero by hand?
    #ifdef SUBZERO
    equate_constant( &Input -> Data.y[shift] , 0.0 ,
		     Input -> Data.y[shift].NSAMPLES ,
		     Input -> Data.y[shift].restype ) ;
    #endif
    
    // compute the partial sum up to whatever cutoff
    for( j = shift+1 ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      
      add( &Input -> Data.y[j] , Input -> Data.y[j-1] ) ;
    }
    
    shift = j ;
  }

  average_equivalent( Input ) ;

  // compute the derivative of the partial sum
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim; i++ ) {

    struct resampled der = init_dist( NULL ,
				      Input -> Data.x[shift].NSAMPLES ,
				      Input -> Data.x[shift].restype ) ;

    for( j = shift+1 ; j < shift + Input -> Data.Ndata[i] - 1 ; j++ ) {

	equate( &der , Input -> Data.y[j+1] ) ;
	subtract( &der , Input -> Data.y[j-1] ) ;
	mult_constant( &der , 1.0/( Input -> Data.x[j+1].avg - Input -> Data.x[j-1].avg ) ) ;
	fprintf( stdout , "[DER] %e %e %e \n" , Input -> Data.x[j].avg , der.avg , der.err ) ;

	if( ( der.avg - der.err ) < 0.0 ) {
	  fprintf( stdout , "CANDIDATE %e %e %e \n" , Input -> Data.x[j].avg ,
		   Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
	  break ;
	}
    }

    free( der.resampled ) ;
    
    shift = j ;
  }

  // multiply
  shift = 0 ;
  const double fac = 1.0 / ( 10*10*10*10. ) ;
  for( i = 0 ; i < Input -> Data.Nsim; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift = j ;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
	
  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}

// fit the slab method evaluation of Qsusc
int
fit_Qslab( struct input_params *Input )
{
  size_t i , j , shift = 0 ;

  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    //const double fac = Input -> Traj[i].Dimensions[3] / ( 10 * 10 * 10 * 10. ) ;

    const double fac =
      //pow( 1.438525 , 4 ) *
      //pow( 1.622301 , 4 ) *
      //pow( 1.808081 , 4 ) *
      //pow( 1.988718 , 4 ) *
      //pow( 2.172728 , 4 ) *
      //pow( 2.356209 , 4 ) *
      //pow( 2.662277 , 4 ) *
      // pow( 2.919241 , 4 ) *
      // pow( 3.799364 , 4 ) *
      //pow( 5.731501 , 4 ) *
      //pow( 7.934524 , 4 ) *
      /*
      Input -> Traj[i].Dimensions[3] / (double)( Input -> Traj[i].Dimensions[0] *
					 Input -> Traj[i].Dimensions[1] *
					 Input -> Traj[i].Dimensions[2] *
					 Input -> Traj[i].Dimensions[3] ) ;
      */
      //1 ;
      Input -> Traj[i].Dimensions[3] ;
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      add_constant( &Input -> Data.x[j] , 1.0 ) ;
      mult_constant( &Input -> Data.y[j] , fac ) ;
    }
    shift=j;
  }

  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    struct resampled res = init_dist( NULL ,
				      Input -> Data.y[shift].NSAMPLES ,
				      Input -> Data.y[shift].restype ) ;
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      if( j == (shift + Input -> Data.Ndata[i] - 1) ) continue ;
      equate( &res , Input -> Data.y[j+1] ) ;
      subtract( &res , Input -> Data.y[j] ) ;
      printf( "%e %e %e \n" , Input -> Data.x[j].avg ,
	      res.avg , res.err ) ;
    }
    shift=j;
  }

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  if( Input -> Fit.Fitdef == QSLAB ) {
    mult_constant( &fit[2] , Input -> Traj[0].Dimensions[3] / 10. ) ;
    
    printf( "T0 * meta :: %e %e \n" , fit[2].avg , fit[2].err ) ;
  }

  // write out a flat file
  char str[256] ;
  sprintf( str , "Qslab_L%zu.flat" , Input -> Traj[0].Dimensions[0] ) ;
  /*
  const double a2 = 100. / ( Input -> Traj[0].Dimensions[0] *
			     Input -> Traj[0].Dimensions[0] ) ;
  */
  FILE *file = fopen( str , "w" ) ;
  //fprintf( file , "%zu\n" , fit[1].restype ) ;
  //fprintf( file , "1\n" ) ;
  fprintf( file , "%zu\n" , fit[1].NSAMPLES ) ;
  for( j = 0 ; j < Input -> Data.y[0].NSAMPLES ; j++ ) {
    fprintf( file , "%1.15e %1.15e\n" , -0.900 , fit[1].resampled[j] ) ;
  }
  fprintf( file , "AVG %1.15e %1.15e\n" , -0.900 , fit[1].avg ) ;
  fclose( file ) ;

  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}

// clean these up
#ifdef SUBZERO
  #undef SUBZERO
#endif
