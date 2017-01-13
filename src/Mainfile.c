/**
   @file Mainfile.c
   @brief Mainly a mainfile
 */
#include "gens.h"

#include "fake.h"
#include "fit_and_plot.h"
#include "init.h"
#include "pmap.h"
#include "read_inputs.h"

int
main( const int argc , const char *argv[] )
{
  size_t i ;
  
  struct input_params Input ;
  if( read_inputs( &Input , "infile" ) == FAILURE ) {
    goto free_failure ;
  }

  // data structure
  Input.Data.Ndata = malloc( Input.Data.Nsim * sizeof( size_t ) ) ;
  Input.Data.Ntot = 0 ;
  for( i = 0 ; i < Input.Data.Nsim ; i++ ) {
    Input.Data.Ndata[i] = 25 ;
    Input.Data.Ntot += Input.Data.Ndata[i] ;
  }
  Input.Data.LT = 25 ;
  
  // need to set this after data has been read ...
  Input.Fit.map = parammap( Input.Data , Input.Fit ) ;

  if( generate_fake_data( &Input.Data , Input.Fit , 0.001 , 0.001 ) == FAILURE ) {
    goto free_failure ;
  }

  if( fit_and_plot( Input ) == FAILURE ) {
    goto free_failure ;
  }

 free_failure :

  // free the structs
  free_Data( &Input.Data , Input.Fit ) ;
  free_Fit( &Input.Fit , Input.Data ) ;
  
  free_inputs( &Input ) ;
  
  return SUCCESS ;
}
