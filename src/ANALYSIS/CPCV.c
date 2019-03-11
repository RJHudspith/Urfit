/**
   @file CPCV.c
   @brief test fitting the CPCV correlators
 */
#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"

int
CVCP_analysis( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  const size_t Ntot_prev = Input -> Data.Ntot ;
  const size_t Nsim_prev = Input -> Data.Nsim ;

  Input -> Data.Ntot /= 2 ;
  Input -> Data.Nsim /= 2 ;

  printf( "What %zu %zu \n" , Input -> Data.Ntot , Input -> Data.Nsim ) ;
  
  // multiply correlator 1 and 2
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      printf( "%e %zu %zu \n" , Input -> Data.y[j].avg , j , Input -> Data.Ntot+j ) ;
      mult( &Input -> Data.y[j] , Input -> Data.y[Input -> Data.Ntot+j] ) ;
      printf( "%e \n" , Input -> Data.y[j].avg ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  
  printf( "Effmass\n" ) ;
  
  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ASINH_EFFMASS ) ;

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  Input -> Data.Ntot = Ntot_prev ;
  Input -> Data.Nsim = Nsim_prev ;
  
  return SUCCESS ;
}
