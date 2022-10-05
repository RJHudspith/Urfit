#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h" // free_fitparams
#include "nrqcd_exp2.h" // set_psq_nrqcd2()
#include "write_flat.h"

#include "resampled_ops.h"

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

  // write out the mass fits
  FILE *file = fopen( "massfits.dat" , "w" ) ;
  struct resampled temp = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;
  struct resampled temp2 = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // p4amp version
    equate( &temp , fit[1] ) ;
    equate( &temp2 , fit[2] ) ;
    raise( &temp2 , -1 ) ;
    mult_constant( &temp2 , p2[i]/2. ) ;
    add( &temp , temp2 ) ;

    fprintf( file , "%f %f\n" , Input->Traj[i].Fit_Low , temp.err_hi ) ;
    fprintf( file , "%f %f\n\n" , Input->Traj[i].Fit_High , temp.err_hi ) ;
    fprintf( file , "%f %f\n" , Input->Traj[i].Fit_Low , temp.avg ) ;
    fprintf( file , "%f %f\n\n" , Input->Traj[i].Fit_High , temp.avg ) ;
    fprintf( file , "%f %f\n" , Input->Traj[i].Fit_Low , temp.err_lo ) ;
    fprintf( file , "%f %f\n\n" , Input->Traj[i].Fit_High , temp.err_lo ) ;    
  }

  fclose( file ) ;

  // write out a flat file of the kinetic mass
  char str[ 256 ] ;
  const double bare_mass = 1.89 ;
  if( Input -> Traj[0].Gs == Vi ) {
    sprintf( str , "Upsilon_%g.flat" , bare_mass ) ;
  } else {
    sprintf( str , "Etab_%g.flat" , bare_mass ) ;
  }

  fprintf( stdout , "Writing mass to %s \n" , str ) ;
  struct resampled mx = init_dist( NULL , fit[2].NSAMPLES , fit[2].restype ) ;
  write_flat_dist( &fit[2] , &mx , 1 , str ) ;

  write_flat_dist( &fit[1] , &mx , 1 , "Mass.flat" ) ;
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
