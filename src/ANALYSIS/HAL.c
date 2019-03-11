#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"

int
HAL_analysis( struct input_params *Input )
{
  // compute the effective mass
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

  size_t i ;
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;
  
  // do the fit
  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    fprintf( stdout , "PHI %e %e %e\n" ,
	     Input -> Traj[i].mom[0]*Input -> Traj[i].mom[0]+
	     Input -> Traj[i].mom[1]*Input -> Traj[i].mom[1]+
	     Input -> Traj[i].mom[2]*Input -> Traj[i].mom[2] ,
	     fit[ 1+i*2 ].avg , fit[ 1+i*2 ].err ) ;
    fprintf( stdout , "MASS %e %e %e\n" ,
	     Input -> Traj[i].mom[0]*Input -> Traj[i].mom[0]+
	     Input -> Traj[i].mom[1]*Input -> Traj[i].mom[1]+
	     Input -> Traj[i].mom[2]*Input -> Traj[i].mom[2] ,
	     fit[ 2+i*2 ].avg , fit[ 2+i*2 ].err ) ;
  }
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
