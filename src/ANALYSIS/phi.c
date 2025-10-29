/**
   @file adler.c
   @brief adler function analysis
 */
#include "gens.h"

#include "stats.h"
#include "resampled_ops.h"
#include "write_flat.h"

static void
compute_phi2( struct resampled rt0 ,
	      struct resampled mpi )
{
  // code assumes sqrt{t_0} is first set to phi
  struct resampled tmp = init_dist( &rt0 ,
				    rt0.NSAMPLES ,
				    rt0.restype ) ;

  printf( "\n t0 %e %e | mpi %e %e\n" , rt0.avg , rt0.err ,
	  mpi.avg , mpi.err ) ;

  mult( &tmp , mpi ) ;
  raise( &tmp , 2 ) ;
  mult_constant( &tmp , 8 ) ;
  printf( "\nPhi2 8.t_0.m^2 :: %f %f\n" , tmp.avg , tmp.err ) ;
  write_flat_dist( &tmp , &rt0 , 1 , "Phi2.flat" ) ;
  free( tmp.resampled ) ;
}

int
phi_analysis( struct input_params *Input )
{
  compute_phi2( Input -> Data.y[0] , Input -> Data.y[1] ) ; 
  
  return SUCCESS ;
}
