/**
   @file adler.c
   @brief adler function analysis
 */
#include "gens.h"

int
write_distributions( struct input_params *Input )
{
  size_t shift = 0 ;
  fprintf( stdout , "\n" ) ;
  for( size_t i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( size_t j = 0 ; j < Input -> Data.Ndata[i] ; j++ ) {
      fprintf( stdout , "%e %e %e\n" ,
	       Input -> Data.x[shift+j].avg ,
	       Input -> Data.y[shift+j].avg ,
	       Input -> Data.y[shift+j].err ) ;
    }
    fprintf( stdout , "\n" ) ;
    shift += Input -> Data.Ndata[i] ;
  }
  
  write_flat_file( *Input , "Resample" ) ;
  return SUCCESS ;
}
