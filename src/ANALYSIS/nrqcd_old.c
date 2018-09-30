/**
   @file nrqcd_old.c
   @brief old-style nrqcd analysis
 */
#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h" // free_fitparams

// just a linear fit
int
nrqcd_old_analysis( struct input_params *Input )
{
  size_t i ;

  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;
  
  // free effective mass
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // do a bunch of single-exponential fits
  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;
  
  for( i = 1 ; i < Input -> Data.Nsim ; i++ ) {

    // subtract the mass at E(0)
    subtract( &fit[2*i+1] , fit[1] ) ;
  }
  subtract( &fit[1] , fit[1] ) ;


  // write a flat file
  char str[ 256 ] ;
  sprintf( str , "masses_ups_psq_%g_%g.flat" , Input -> Traj[0].Fit_Low ,
	   Input -> Traj[0].Fit_High ) ;
  FILE *outfile = fopen( str , "w" ) ;

  fprintf( outfile , "%d\n" , fit[0].restype ) ;
  fprintf( outfile , "%zu\n" , Input -> Data.Nsim ) ;

  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // compute momentum
    size_t mu ;
    double p2 = 0 ;
    for( mu = 0 ; mu < 3 ; mu++ ) {
      const double pref = ( Input -> Traj[i].mom[mu] * 2. * M_PI /
			    Input -> Traj[i].Dimensions[mu] ) ;
      p2 += pref * pref ;
    }
    fprintf( outfile , "%zu\n" , fit[0].NSAMPLES ) ;

    size_t k ;
    for( k = 0 ; k < fit[0].NSAMPLES ; k++ ) {
      fprintf( outfile , "%1.12f %1.12f\n" ,
	       p2 , fit[2*i+1].resampled[k] ) ;
    }
    fprintf( outfile , "AVG %1.12f %1.12f\n" , p2 , fit[2*i+1].avg ) ;
  }

  fclose( outfile ) ;
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  
  return SUCCESS ;
}
