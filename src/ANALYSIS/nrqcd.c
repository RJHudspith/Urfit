#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h" // free_fitparams
#include "nrqcd_exp2.h" // set_psq_nrqcd2()

// just a linear fit
int
nrqcd_analysis( struct input_params *Input )
{
  size_t i ;

  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;
  
  // free effective mass
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // init p^2
  double p2[ Input -> Data.Nsim ] ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t mu ;
    p2[i] = 0 ;
    for( mu = 0 ; mu < 3 ; mu++ ) {
      const double pref = ( Input -> Traj[i].mom[mu] * 2. * M_PI /
			    Input -> Traj[i].Dimensions[mu] ) ;
      p2[i] += pref * pref ;
    }
  }

  // set the psq array
  set_psq_nrqcd2( p2 , Input -> Data.Nsim ) ;
  
  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  // write out a flat file of the kinetic mass
  char str[ 256 ] ;
  const double bare_mass = 9.0 ;
  if( Input -> Traj[0].Gs == Vi ) {
    sprintf( str , "Upsilon_%g.flat" , bare_mass ) ;
  } else {
    sprintf( str , "Etab_%g.flat" , bare_mass ) ;
  }

  fprintf( stdout , "Writing mass to %s \n" , str ) ;

  FILE *outfile = fopen( str , "w" ) ;
  fprintf( outfile , "%d\n" , fit[2].restype ) ;
  fprintf( outfile , "1\n" ) ;
  fprintf( outfile , "%zu\n" , fit[2].NSAMPLES ) ;
  
  for( i = 0 ; i < fit[2].NSAMPLES ; i++ ) {
    fprintf( outfile , "%1.15e %1.15e\n" , bare_mass , fit[2].resampled[i] ) ;
  }
  fprintf( outfile , "AVG %1.15e %1.15e\n" , bare_mass , fit[2].avg ) ;
  
  fclose( outfile ) ;
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
