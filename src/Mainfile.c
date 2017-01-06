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
#include "stats.h"

// minimizers
#include "GA.h"
#include "CG.h"
#include "SD.h"
#include "LM.h"

int main( void )
{
  size_t i ;

  // data structure
  struct data_info Data ;
  Data.Nsim = 3 ;
  Data.Nboots = 500 ;
  Data.Ndata = malloc( Data.Nsim * sizeof( size_t ) ) ;
  Data.Ntot = 0 ;
  for( i = 0 ; i < Data.Nsim ; i++ ) {
    Data.Ndata[i] = 35 ;
    Data.Ntot += Data.Ndata[i] ;
  }
  Data.LT = 25 ;
  Data.Cov.Divided_Covariance = false ;
  Data.Cov.Column_Balanced = false ;
  Data.Cov.Eigenvalue_Tol = 1E-8 ;

  // fit structure
  struct fit_info Fit ;
  Fit.Fitdef = PADE ;
  Fit.Corrfit = UNCORRELATED ;
  Fit.N = 2 ;
  Fit.M = 1 ;
  Fit.Nparam = get_Nparam( Fit ) ;
  Fit.Sims = malloc( Fit.Nparam * sizeof( bool ) ) ;
  for( i = 0 ; i < Fit.Nparam ; i++ ) {
    Fit.Sims[i] = false ; i&1 ? false : true ;
  }
  Fit.map = parammap( Data , Fit ) ;
  Fit.Minimize = ga_iter ;
  Fit.Tol = 1E-7 ;

  struct resampled *fitparams = NULL ;

  if( generate_fake_data( &Data , Fit , 0.001 , 0.001 ) == FAILURE ) {
    goto free_failure ;
  }

  if( inverse_correlation( &Data , Fit ) == FAILURE ) {
    goto free_failure ;
  }

#ifdef VERBOSE
  write_corrmatrix( (const double**)Data.Cov.W , Data.Ntot ) ;
#endif
  
  if( ( fitparams = perform_bootfit( Data , Fit ) ) == NULL ) {
    goto free_failure ;
  }

  make_graph( fitparams , Data , Fit , "test.agr" , "x" , "y" ) ;

 free_failure :

  // free the fit
  if( fitparams != NULL ) {
    for( i = 0 ; i < Fit.Nparam ; i++ ) {
      if( fitparams[i].resampled != NULL ) {
	free( fitparams[i].resampled ) ;
      }
    }
    free( fitparams ) ;
  }

  free_Data( &Data ) ;
  free_Fit( &Fit , Data ) ;
  
  return SUCCESS ;
}
