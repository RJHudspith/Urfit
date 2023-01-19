/**
   Spin-Orbit calculation
 */
#include "gens.h"

#include "stats.h"
#include "resampled_ops.h"

// does (3*mV+mP )/4.
void
SO_split0( const struct resampled mP ,
	   const struct resampled mV ,
	   const double ainv )
{
  struct resampled temp = init_dist( NULL , mV.NSAMPLES , mV.restype ) ;
  for( size_t i = 0 ; i < mV.NSAMPLES ; i++ ) {
    temp.resampled[i] = ( 3*mV.resampled[i] + mP.resampled[i] )/4. ;
  }
  temp.avg = ( 3*mV.avg + mP.avg )/4. ;
  compute_err( &temp ) ;
  fprintf( stdout , "Spin-avg %e %e\n" ,
	   temp.avg*ainv , temp.err*ainv ) ;
  free( temp.resampled ) ;
}

// does mV-mP
void
hyperfine1( const struct resampled mP ,
	    const struct resampled mV ,
	    const double ainv )
{
  struct resampled temp = init_dist( NULL , mV.NSAMPLES , mV.restype ) ;
  for( size_t i = 0 ; i < mV.NSAMPLES ; i++ ) {
    temp.resampled[i] = ( mV.resampled[i] - mP.resampled[i] ) ;
  }
  temp.avg = ( mV.avg - mP.avg ) ;
  compute_err( &temp ) ;
  fprintf( stdout , "Hyperfine %e %e\n" ,
	   temp.avg*ainv , temp.err*ainv ) ;
  free( temp.resampled ) ;
}

// does (5*mT-3*mA-2*mI)/9
void
SO_split1( const struct resampled mI ,
	   const struct resampled mA ,
	   const struct resampled mT ,
	   const double ainv )
{
  struct resampled temp = init_dist( NULL , mI.NSAMPLES , mI.restype ) ;
  for( size_t i = 0 ; i < mI.NSAMPLES ; i++ ) {
    temp.resampled[i] = ( +5*mT.resampled[i]
			  -3*mA.resampled[i]
			  -2*mI.resampled[i] ) /9. ;
  }
  temp.avg = ( 5*mT.avg - 3*mA.avg - 2*mI.avg )/9. ;
  compute_err( &temp ) ;
  fprintf( stdout , "Spin-orbit %e %e\n" ,
	   temp.avg*ainv , temp.err*ainv ) ;
  free( temp.resampled ) ;
}

// does (3*mA-mT-2*mI)/9
void
SO_split2( const struct resampled mI ,
	   const struct resampled mA ,
	   const struct resampled mT ,
	   const double ainv )
{
  struct resampled temp = init_dist( NULL , mI.NSAMPLES , mI.restype ) ;
  for( size_t i = 0 ; i < mI.NSAMPLES ; i++ ) {
    temp.resampled[i] = ( -mT.resampled[i]
			  +3*mA.resampled[i]
			  -2*mI.resampled[i] ) /9. ;
  }
  temp.avg = ( -mT.avg + 3*mA.avg - 2*mI.avg )/9. ;
  compute_err( &temp ) ;
  fprintf( stdout , "tensor split %e %e\n" ,
	   temp.avg*ainv , temp.err*ainv ) ;
  free( temp.resampled ) ;
}

// does (5*mT+3*mT+2*mI)/9
void
SO_split3( const struct resampled mI ,
	   const struct resampled mA ,
	   const struct resampled mT ,
	   const double ainv )
{
  struct resampled temp = init_dist( NULL , mI.NSAMPLES , mI.restype ) ;
  for( size_t i = 0 ; i < mI.NSAMPLES ; i++ ) {
    temp.resampled[i] = ( +5*mT.resampled[i]
			  +3*mA.resampled[i]
			  +2*mI.resampled[i] ) /9. ;
  }
  temp.avg = ( 5*mT.avg + 3*mA.avg + 2*mI.avg )/9. ;
  compute_err( &temp ) ;
  fprintf( stdout , "1P spinavg %e %e\n" ,
	   temp.avg*ainv , temp.err*ainv ) ;
  free( temp.resampled ) ;
}

int
spin_orbit( struct input_params *Input )
{
  if( Input -> Data.Nsim != 5 ) {
    fprintf( stdout , "[SpinOrbit] requires 5 arguments: eta_c, j_psi, scalar, axial, tensor\n" ) ;    
    return FAILURE ;
  }

  const double hc = 0.1973269804 ; 
  //const double ainv = (hc/0.09929) ;
  //const double ainv = (hc/0.08636) ;
  //const double ainv = (hc/0.07634) ;
  //const double ainv = (hc/0.06426) ;
  const double ainv = (hc/0.04981) ;

  // print the results
  fprintf(stdout , "\nMASSES\n" ) ;
  fprintf( stdout , "Metac %e %e\n" , Input -> Data.y[0].avg*ainv , Input -> Data.y[0].err*ainv ) ;
  fprintf( stdout , "MJPsi %e %e\n" , Input -> Data.y[1].avg*ainv , Input -> Data.y[1].err*ainv ) ;
  fprintf( stdout , "Mc0 %e %e\n" , Input -> Data.y[2].avg*ainv , Input -> Data.y[2].err*ainv ) ;
  fprintf( stdout , "Mc1 %e %e\n" , Input -> Data.y[3].avg*ainv , Input -> Data.y[3].err*ainv ) ;
  fprintf( stdout , "Mc2 %e %e\n" , Input -> Data.y[4].avg*ainv , Input -> Data.y[4].err*ainv ) ;

  // spin averages
  fprintf(stdout , "\nSPINAVG\n" ) ;
  SO_split0( Input -> Data.y[0] , Input -> Data.y[1] , ainv ) ;
  SO_split3( Input -> Data.y[2] , Input -> Data.y[3] , Input -> Data.y[4] , ainv ) ;

  // spin splittings
  fprintf(stdout , "\nSPINSPLIT\n" ) ;
  hyperfine1( Input -> Data.y[0] , Input -> Data.y[1] , ainv ) ;
  SO_split1( Input -> Data.y[2] , Input -> Data.y[3] , Input -> Data.y[4] , ainv ) ;
  SO_split2( Input -> Data.y[2] , Input -> Data.y[3] , Input -> Data.y[4] , ainv ) ;
  
  return 0 ;
}
