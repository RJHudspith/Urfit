#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h" // free_fitparams

#include "fsol2.h" // set_psq_sol()
//#include "fsol.h" // set_psq_sol()

#include "resampled_ops.h"

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
      double pref = ( Input -> Traj[i].mom[mu] * 2. * M_PI /
		      Input -> Traj[i].Dimensions[mu] ) ;
      //pref = 2*sin(pref*0.5) ;
      p2[i] += pref * pref ;
    }
  }

  // set the psq array
  set_psq_sol2( p2 , Input -> Data.Nsim ) ;
  
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


  // write a flat file
  equate_constant( &temp , pow(1./Input->Traj[0].Dimensions[2],2) ,
		   temp.NSAMPLES , temp.restype ) ;  
  mult_constant( &fit[1] , Input->Traj[0].Dimensions[2] ) ;
  write_flat_dist( &fit[1] , &temp , 1 , "m0.flat" ) ;

  // write flat disperion file
  equate_constant( &temp , pow(1./Input->Traj[0].Dimensions[2],2) ,
		   temp.NSAMPLES , temp.restype ) ;  
  // mult_constant( &fit[5] , Input->Traj[0].Dimensions[2] ) ;
  write_flat_dist( &fit[5] , &temp , 1 , "csq.flat" ) ;
  
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
