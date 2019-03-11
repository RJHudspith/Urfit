/**
   @file hvp_pade.c
   @brief fit a pade to the hvp
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "fit_chooser.h"
#include "ffunction.h"
#include "init.h"
#include "make_xmgrace.h"
#include "resampled_ops.h"
#include "stats.h"

int
su2_shit( struct input_params *Input )
{
  size_t i = 0 ;
  double chisq ;

  const double xmap[10] = { 0.2606 , 0.2036 , 0.2290 ,
			    0.1827 , 0.1657 , 0.1179 ,
			    0.0793 , 0.0578 , 0.0444 ,
			    0.0353 } ;
  
  size_t j , shift = 0 ;
  // take the log of both sides
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {

    equate_constant( &Input -> Data.x[i] , xmap[i%10] ,
		     Input -> Data.x[i].NSAMPLES , Input -> Data.x[i].restype ) ;
    mult_constant( &Input -> Data.y[i] , 1.3317 ) ;

    
    //raise( &Input -> Data.y[i] , 2 ) ;
    //printf( "%1.12e %1.12e\n" , Input -> Data.y[i].avg , Input -> Data.y[i].err ) ;
    /*
    const double V = 1.0 / ( Input -> Traj[i].Dimensions[0] *
			     Input -> Traj[i].Dimensions[1] *
			     Input -> Traj[i].Dimensions[2] *
			     Input -> Traj[i].Dimensions[3] ) ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      mult_constant( &Input -> Data.y[j] , V ) ;
    }
    shift += Input -> Data.Ndata[i] ;
    */
    //mult_constant( &Input -> Data.y[i] , -1 ) ;
    //add_constant( &Input -> Data.x[i] , 0.0019 ) ;
  }

  for( j = 0 ; j < 10 ; j++ ) {
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      fprintf( stdout , "%e %1.15e %1.15e\n" ,
	       1.3317*1.3317*1.3317/(Input -> Traj[i].Dimensions[0]*
		   Input -> Traj[i].Dimensions[1]*
		   Input -> Traj[i].Dimensions[2]) ,
	       Input -> Data.y[j+10*i].avg ,
	       Input -> Data.y[j+10*i].err ) ;
    }
    printf( "\n" ) ;
  }
  
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  struct resampled res = init_dist( NULL ,
				    fit[0].NSAMPLES ,
				    fit[0].restype ) ;

  /*
  subtract( &res , fit[0] ) ;
  divide( &res , fit[1] ) ;

  printf( "Zero %e %e \n" , res.avg , res.err ) ;

  mult_constant( &res , 0.0 ) ;
  subtract( &res , fit[2] ) ;
  divide( &res , fit[3] ) ;

  printf( "Zero %e %e \n" , res.avg , res.err ) ;
  */

  //fprintf( stdout , "PARAM %e %e\n" , fit[0].avg , fit[0].err ) ;
  //root( &fit[0] ) ;
  //fprintf( stdout , "PARAM %e %e\n" , fit[0].avg , fit[0].err ) ;
  
  return SUCCESS ;
}
