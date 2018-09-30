#include "gens.h"
#include "resampled_ops.h"
#include "fit_and_plot.h"

static const double c4map[12] = { 1,1.23,1.4 ,
				  1,1.23,1.4 ,
				  1,1.23,1.4 ,
				  1,1.23,1.4 } ;

int
c4c7_analysis( struct input_params *Input )
{
  size_t i ;
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    equate_constant( &Input -> Data.x[i] , c4map[i] ,
		     Input -> Data.x[i].NSAMPLES ,
		     Input -> Data.x[i].restype ) ;
    printf( "HERE %f \n" , Input -> Data.x[i].avg ) ;
  }

  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
