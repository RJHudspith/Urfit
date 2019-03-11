#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "momenta.h"
#include "resampled_ops.h"
#include "write_flat.h"

// topological moments of Q^2
int
Qmoments( struct input_params *Input )
{
  size_t j ;

  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    printf( "[QMOM] QMOM_%zu %f %f \n" , j ,
	    Input -> Data.y[j].avg ,
	    Input -> Data.y[j].err ) ;
  }
  
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    if( j!=0 ) {
      divide( &Input -> Data.y[j] , Input -> Data.y[0] ) ;
    }
  }

  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    printf( "[QMOM] QNORM_%zu %f %f \n" , j ,
	    Input -> Data.y[j].avg ,
	    Input -> Data.y[j].err ) ;
  }

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  // write out some files
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    char str[ 256 ] ;
    sprintf( str , "Qmoment.%zu.L%zu.flat" ,
	     j , Input -> Traj[0].Dimensions[0] ) ;
    const double a2 = 10. * 10. /
      ( Input -> Traj[0].Dimensions[0] *
	Input -> Traj[0].Dimensions[0] ) ;
    
    struct resampled x = init_dist( NULL ,
				    Input -> Data.y[j].NSAMPLES ,
				    Input -> Data.y[j].restype ) ;
    equate_constant( &x , a2 , x.NSAMPLES , x.restype ) ;
    write_flat_dist( &Input -> Data.y[j] , &x , 1 , str ) ;

    free( x.resampled ) ;
  }
  
  return SUCCESS ;
}
