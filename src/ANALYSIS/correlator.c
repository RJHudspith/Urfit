/**
   @file correlator.c
   @brief correlator analysis
 */
#include "gens.h"

#include "decays.h"
#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"

//#define FIT_EFFMASS
//#define COMBINE

int
correlator_analysis( struct input_params *Input )
{
  size_t i ;
#ifdef COMBINE
  for( i = 0 ; i < Input -> Data.Ndata[0] ; i++ ) {
    divide_constant( &Input -> Data.y[i+Input -> Data.Ndata[0]] , 12*12*12*32 ) ;
    subtract( &Input -> Data.y[i] , Input -> Data.y[i+Input -> Data.Ndata[0]] ) ;
    equate( &Input -> Data.y[i+Input -> Data.Ndata[0]] , Input -> Data.y[i] ) ;
  }
#endif

  printf( "Effmass\n" ) ;
  
  // compute an effective mass ? TODO
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
    
  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  // compute a decay constant
  if( Input -> Fit.Fitdef == PP_AA_WW ||
      Input -> Fit.Fitdef == PP_AA ||
      Input -> Fit.Fitdef == PPAA ) {
    struct resampled dec = decay( Fit , *Input ) ;

    divide( &Fit[0] , dec ) ;

    raise( &Fit[0] , 2 ) ;
    
    printf( "(M/F)^2 %e %e \n" , Fit[0].avg , Fit[0].err ) ; 
    
    free( dec.resampled ) ;
  }  

  /*
  raise( &Fit[3] , 2 ) ;

  printf( "M^2 :: %f %f \n" , Fit[3].avg , Fit[3].err ) ;
  */
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
