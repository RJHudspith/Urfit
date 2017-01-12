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

  struct data_info Data ;
  struct resampled *fitparams = NULL ;
  
  struct input_params Input ;
  if( read_inputs( &Input , "infile" ) == FAILURE ) {
    goto free_failure ;
  }

  // data structure
  Data.Nsim = Input.Ntraj ;
  Data.Nboots = 500 ;
  Data.Ndata = malloc( Data.Nsim * sizeof( size_t ) ) ;
  Data.Ntot = 0 ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    Data.Ndata[i] = 25 ;
    Data.Ntot += Data.Ndata[i] ;
  }
  Data.LT = 25 ;
  Data.Cov.Divided_Covariance = false ;
  Data.Cov.Column_Balanced = false ;
  Data.Cov.Eigenvalue_Tol = 1E-8 ;
  
  // need to set this after data has been read ...
  Input.Fit.map = parammap( Data , Input.Fit ) ;

  if( generate_fake_data( &Data , Input.Fit , 0.001 , 0.001 ) == FAILURE ) {
    goto free_failure ;
  }

  if( inverse_correlation( &Data , Input.Fit ) == FAILURE ) {
    goto free_failure ;
  }
  
#ifdef VERBOSE
  write_corrmatrix( (const double**)Data.Cov.W , Data.Ntot ) ;
#endif
  
  if( ( fitparams = perform_bootfit( Data , Input.Fit ) ) == NULL ) {
    goto free_failure ;
  }

  make_graph( fitparams , Data , Input.Fit , "test.agr" , "x" , "y" ) ;

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
  free_Data( &Data ) ;
  free_Fit( &Input.Fit , Data ) ;
  free_inputs( &Input ) ;
  
  return SUCCESS ;
}
