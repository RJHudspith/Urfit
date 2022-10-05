#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h" // free_fitparams
#include "fsol.h" // set_psq_sol()

#include "resampled_ops.h"

// just a linear fit
int
sol_analysis( struct input_params *Input )
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
  set_psq_sol( p2 , Input -> Data.Nsim ) ;
  
  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  // write out the mass fits
  FILE *file = fopen( "massfits.dat" , "w" ) ;
  struct resampled temp = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;
  struct resampled temp2 = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    equate( &temp , fit[1] ) ;
    raise( &temp , 2 ) ;
    
    equate( &temp2 , fit[2] ) ;
    mult_constant( &temp2 , p2[i] ) ;

    add( &temp , temp2 ) ;
    root( &temp ) ;

    fprintf( file , "%f %f\n" , Input->Traj[i].Fit_Low , temp.err_hi ) ;
    fprintf( file , "%f %f\n\n" , Input->Traj[i].Fit_High , temp.err_hi ) ;
    fprintf( file , "%f %f\n" , Input->Traj[i].Fit_Low , temp.avg ) ;
    fprintf( file , "%f %f\n\n" , Input->Traj[i].Fit_High , temp.avg ) ;
    fprintf( file , "%f %f\n" , Input->Traj[i].Fit_Low , temp.err_lo ) ;
    fprintf( file , "%f %f\n\n" , Input->Traj[i].Fit_High , temp.err_lo ) ;    
  }

  fclose( file ) ;

  // write out a flat file of the kinetic mass
  char str[ 256 ] , strD[ 256 ] ;
  const double bare_mass = 1.89 ;
  if( Input -> Traj[0].Gs == Vi ) {
    sprintf( str , "Upsilon.flat" , bare_mass ) ;
    sprintf( strD , "Upsilon_disp.flat" , bare_mass ) ;
  } else {
    sprintf( str , "Etab.flat" , bare_mass ) ;
    sprintf( strD , "Etab_disp.flat" , bare_mass ) ;
  }

  fprintf( stdout , "Writing mass to %s \n" , str ) ;
  fprintf( stdout , "Writing disp to %s \n" , strD ) ;

  FILE *outfile  = fopen( str  , "w" ) ;
  FILE *outfileD = fopen( strD , "w" ) ;
  fprintf( outfile , "%d\n" , fit[2].restype ) ;
  fprintf( outfile , "1\n" ) ;
  fprintf( outfile , "%zu\n" , fit[2].NSAMPLES ) ;
  fprintf( outfileD , "%d\n" , fit[2].restype ) ;
  fprintf( outfileD , "1\n" ) ;
  fprintf( outfileD , "%zu\n" , fit[2].NSAMPLES ) ;

  
  for( i = 0 ; i < fit[2].NSAMPLES ; i++ ) {
    fprintf( outfile , "%1.15e %1.15e\n" , bare_mass , fit[1].resampled[i] ) ;
    fprintf( outfileD , "%1.15e %1.15e\n" , bare_mass , fit[2].resampled[i] ) ;
  }
  fprintf( outfile , "AVG %1.15e %1.15e\n" , bare_mass , fit[1].avg ) ;
  fprintf( outfileD , "AVG %1.15e %1.15e\n" , bare_mass , fit[2].avg ) ;
  
  fclose( outfile ) ;
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
