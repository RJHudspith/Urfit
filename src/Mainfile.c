/**
   @file Mainfile.c
   @brief Mainly a mainfile
 */
#include "gens.h"

#include "bootfit.h"
#include "correlation.h"
#include "fake.h"
#include "fit_chooser.h"
#include "init.h"
#include "make_xmgrace.h"
#include "plot_fitfunc.h"
#include "pmap.h"
#include "read_inputs.h"
#include "stats.h"

// minimizers
#include "GA.h"
#include "CG.h"
#include "SD.h"
#include "LM.h"

int
main( const int argc , const char *argv[] )
{
  size_t i ;

  struct resampled *fitparams = NULL ;
  
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

  if( generate_fake_data( &Input.Data , Input.Fit , 0.1 , 0.1 ) == FAILURE ) {
    goto free_failure ;
  }

  if( inverse_correlation( &Input.Data , Input.Fit ) == FAILURE ) {
    goto free_failure ;
  }
  
  //#ifdef VERBOSE
  write_corrmatrix( (const double**)Input.Data.Cov.W ,
		    Input.Data.Ntot , Input.Fit.Corrfit ) ;
  //#endif
  
  if( ( fitparams = perform_bootfit( Input.Data , Input.Fit ) ) == NULL ) {
    goto free_failure ;
  }

  // and make a graph!
  make_graph( fitparams , Input.Data , Input.Fit , Input.Graph ) ;

 free_failure :

  // free the fit
  if( fitparams != NULL ) {
    for( i = 0 ; i < Input.Fit.Nparam ; i++ ) {
      if( fitparams[i].resampled != NULL ) {
	free( fitparams[i].resampled ) ;
      }
    }
    free( fitparams ) ;
  }

  // free the structs
  free_Data( &Input.Data , Input.Fit ) ;
  free_Fit( &Input.Fit , Input.Data ) ;
  
  free_inputs( &Input ) ;
  
  return SUCCESS ;
}
