#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"

int
binding_corr_analysis( struct input_params *Input )
{
  size_t i , j ;
  const size_t N = Input -> Data.Ndata[0] ;

  for( i = 0 ; i < Input -> Data.Nsim/3 ; i++ ) {
    for( j = 0 ; j < N ; j++ ) {
      mult( &Input -> Data.y[ j+(3*i)*N ] ,
	    Input -> Data.y[ j+(1+3*i)*N ] ) ;
      divide( &Input -> Data.y[ j+(2+3*i)*N ] ,
	      Input -> Data.y[ j+(3*i)*N ] ) ;
    }
  }
  
  fprintf( stdout , "Effmass\n" ) ;
  
  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}